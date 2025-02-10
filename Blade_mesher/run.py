from IPython import get_ipython
get_ipython().run_line_magic('clear', '')

from timeit import default_timer as timer
start_time = timer()  # Start the timer

import numpy as np, os, sys 

# Get the absolute path to the directory and add it to the sys path for relative imports
package_dir = os.path.abspath(os.path.dirname(__file__))
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
n_blades = 1        # Number of blades

remove_TE = True            # Should we remove TE ??
close_TE = True             # Should we close the TE gap ??

######### .pmt Inputs #########
surf_type = 'propeller'     # Surface type (wing/propeller)
rotation_center = [0, 0, 0]
rotation_axis = [0, -1, 0]

############################ Panel Parameters #################################
### Airfoil Section ###
num_points = 25                         ## Upper and Lower surf separately!!
dist_airfoil = 'cosine_LE'

### Spanwise Cutting Planes ###
spanwise_panel_num= 55
z_min, z_max = 24, 152                ## in mm

dist_spanwise = 'cosine_TIP'
r_R = 0.83                               ## Span location to start the cosine_TIP 
"""###############################################################################"""
if close_TE is True:
    TE_property = ''
else:
    TE_property='-OPEN_TE'    
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
# Generate mesh #
mesh, connectivity, mesh_DUST, coordinates_DUST, connectivity_DUST = generate_mesh(ALL_sections, n_blades, close_TE)

import meshio
meshio.write(f'{output_dir}/12x6_mesh___span_{dist_spanwise}-sec_{dist_airfoil}{TE_property}.vtk', mesh, file_format='vtk')

"""###############  Export the point coordinates and their connectivity information (for Basic Mesh in DUST)   ##################"""
DUST_dir = f'{output_dir}/DUST_input___{len(ALL_sections)}_n_span-{ALL_sections[0].shape[0]}_n_sec{TE_property}'
if os.path.isdir(DUST_dir) is False:
    os.mkdir(f'{DUST_dir}')

meshio.write(f'{DUST_dir}/mesh_DUST.vtk', mesh_DUST, file_format='vtk')

with open(f'{DUST_dir}/rr.dat', 'w') as file:
    np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

with open(f'{DUST_dir}/ee.dat', 'w') as file:
    np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')
"""##############################################################################################################################"""    

# ## Export all cross sections (Optional) ##
# cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")

######### Write the all coordinates to a .pmt file to run with the code #######

with open(f'{output_dir}/12x6_Blade___span_{dist_spanwise}-sec_{dist_airfoil}{TE_property}.pmt', 'w') as file:
    file.write('######## Panel parameters ########\n')
    file.write(f'type={surf_type}\n')
    file.write(f'n_blades={n_blades}\n\n')
    file.write(f'rotation_center={rotation_center}\n')
    file.write(f'rotation_axis={rotation_axis}\n')
    
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

###############################################################################################################################
end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds")