from IPython import get_ipython
get_ipython().run_line_magic('clear', '')

from timeit import default_timer as timer
start_time = timer()  # Start the timer

import numpy as np, os, sys 
import cadquery as cq
import meshio

# Get the absolute path to the parent directory and add it to the sys path for relative imports
package_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if package_dir not in sys.path:
    sys.path.append(package_dir)

from Distributions.distributions import spanwise_disribution as spanwise_planes
from Point_coordinates import point_generation as points
from Mesher.mesh import generate_mesh

"""################################ INPUTS #######################################"""
# File paths and the STEP file
working_dir = '../Test_Cases/12x6 ClarkY'
stp_file = f'{working_dir}/12x6_ClarkY-1_Blade_LE_at_the_mid_face_little_modification.stp'

############################ Geometry Parameters ##############################
number_of_blades = 2        # Number of blades

remove_TE = True            # Should we remove TE ??
close_TE = True             # Should we close the TE gap ??

######### .pmt Inputs #########
surf_type = 'propeller'     # Surface type (wing/propeller)
rotation_center = [0, 0, 0]
rotation_axis = [0, -1, 0]

############################ Panel Parameters #################################
### Airfoil Section ###
num_points = 19                         ## Upper and Lower surf separately!!
dist_airfoil = 'cosine_LE'

### Spanwise Cutting Planes ###
spanwise_panel_num= 52
z_min, z_max = 26, 152                ## in mm

dist_spanwise = 'cosine_TIP'
r_R = 0.8                               ## Span location to start the cosine_TIP 

"""###############################################################################"""

if close_TE is True:
    TE_property = '-closed_TE'
else:
    TE_property=''
################################

output_dir = f'{working_dir}/output'
# Create output directory if there is no #
if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

### Create spanwise planes ###
z_planes = spanwise_planes(z_min, z_max, spanwise_panel_num, dist_spanwise, r_R)

###### GENERATE POINTS ######
ALL_sections, all_sections_compound = points.get_coords(stp_file, num_points, dist_airfoil, z_planes, remove_TE, close_TE)
    
###### GENERATE MESH AND EXPORT #######
mesh = generate_mesh(ALL_sections, number_of_blades)
meshio.write(f'{output_dir}/12x6_mesh___span_{dist_spanwise}-sec_{dist_airfoil}{TE_property}.cgns', mesh, file_format='cgns')


# ## Export all cross sections (Optional) ##
# cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")

######### Write the all coordinates to a .pmt file to run with the code #######

with open(f'{output_dir}/12x6_Blade___span_{dist_spanwise}-sec_{dist_airfoil}{TE_property}.pmt', 'w') as file:
    file.write('######## Panel parameters ########\n')
    file.write('type=' + surf_type + '\n')
    file.write('n_blades=' + str(number_of_blades) + '\n\n')
    file.write('rotation_center=' + str(rotation_center) + '\n')
    file.write('rotation_axis=' + str(rotation_axis) + '\n')
    
    file.write('n_span_all=' + str(len(ALL_sections)) + '\n')
    file.write('n_points=' + str(len(ALL_sections[0][:,0])) + '\n')
    file.write('######## End of parameters ########\n')
    
    for i in range(len(ALL_sections)):
        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')
        
##### Write coordinates to txt to check with the paraview ######
with open(f'{output_dir}/Blade_points_check.txt', 'w') as file:
    file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
    for i in range(len(ALL_sections)):
        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

"""
# Propeller coordinates
with open(f'{output_dir}/propeller_points.txt', 'w') as file:
    file.write('X(mm)\tY(mm)\tZ(mm)\n')
    np.savetxt(file, mesh.points, delimiter='\t', fmt=['%.9f','%.9f','%.9f'], comments='')
"""

###############################################################################################################################
end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds")