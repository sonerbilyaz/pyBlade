from IPython import get_ipython
get_ipython().run_line_magic('clear', '')

from timeit import default_timer as timer
start_time = timer()  # Start the timer

from Distributions.distributions import airfoil_distribution as airfoil_points
from Distributions.distributions import spanwise_disribution as spanwise_planes

import Geometry_operations.modify as modify
import Geometry_operations.extract_info as identify

import cadquery as cq
import numpy as np, os 

"""################################ INPUTS #######################################"""
# File paths and the STEP file
working_dir = 'Test_Cases/12x6 ClarkY'
stp_file = f'{working_dir}/12x6_ClarkY-1_Blade_LE_at_the_mid_face_little_modification.stp'

############################ Geometry Parameters ##############################
surf_type = 'propeller'     # Surface type (wing/propeller)
number_of_blades = 2        # Number of blades

remove_TE = True            # Should we remove TE ??
close_TE = True             # Should we close the TE gap ??

rotation_center = [0, 0, 0]
rotation_axis = [0, -1, 0]

############################ Panel Parameters ##################################
""" Airfoil Section """
# Number of points to generate (upper and lower surf separately)
num_points = 19
dist_airfoil = 'cosine_LE'

# TE_inflation = 0.1  ## r/R ratio at which TE detection method will switch (due to the local pitch change from root to tip) based on the x-coordinate comparison

""" Spanwise Cutting Planes """
spanwise_panel_num= 50
z_min, z_max = 24, 152              ## in mm
dist_spanwise = 'cosine_TIP'
r_R = 0.75                           ## Start span location of the cosine_TIP 
  
"""###############################################################################"""

output_dir = f'{working_dir}/output'
# Create output directory if there is no #
if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

### Create parametric points for airfoil distribution (FOR UPPER&LOWER SURFACES SEPERATELY) ###
parametric_points, Node_ID_one_surf = airfoil_points(num_points, dist_airfoil)

### Create spanwise planes ###
z_planes = spanwise_planes(z_min, z_max, spanwise_panel_num, dist_spanwise, r_R)

## Get coordinates ##
def get_coords(stp_file, num_points, z_planes, remove_TE):    
    # Load the .stp file
    blade = cq.importers.importStep(stp_file)
    
    ALL_sections = []
    ALL_cross_sections = []
    
    global z
    for z in z_planes:
        # Create a section of the blade using the defined z-coord
        cross_section=blade.section(height=z)
        # Get all edges from the cross section
        edges = cross_section.edges()
        # Access the edge objects from the Workplane
        edge_objects = edges.objects
        
        # Should we remove TE ??
        if remove_TE is True:
            filtered_edges, TE_upper, TE_lower = modify.remove_TE(edge_objects)
        
        if remove_TE is False:
            _, TE_upper, TE_lower  = modify.remove_TE(edge_objects)
            filtered_edges = edge_objects
        
        # Identify the edges
        upper_surface, lower_surface, LE_coords = identify.extract_edges(filtered_edges, TE_upper, TE_lower)
        
        ## Assign the parametric points between 0 and 1 to the upper and lower surfaces
        parametric_points_up = parametric_points 
        Node_ID_up = Node_ID_one_surf
        
        parametric_points_low = parametric_points
        Node_ID_low = np.linspace(len(Node_ID_up)+1, 2*len(Node_ID_up)-1, len(Node_ID_up)-1)
        
        # Check it starts generating points from the starting vertex (LE)
        if np.array_equal(np.round(upper_surface.positionAt(0).toTuple(),decimals=2).reshape(1,3), np.round(LE_coords,decimals=2)) is False and dist_airfoil == 'cosine_LE':
            # If the orientation is wrong, reverse it to correct ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
            parametric_points_up = 1 - parametric_points_up
            # 
            # print(f'upper surface parametric is reversed at z={z}mm. LE = {upper_surface.positionAt(parametric_points_up[0]).toTuple()}')
            
        if np.array_equal(np.round(lower_surface.positionAt(0).toTuple(),decimals=2).reshape(1,3), np.round(LE_coords,decimals=2)) is False and dist_airfoil == 'cosine_LE':
            # If the orientation is wrong, reverse it to correct ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
            parametric_points_low = 1- parametric_points_low
            # 
            # print(f'lower surface parametric is reversed at z={z}mm LE = {lower_surface.positionAt(parametric_points_low[0]).toTuple()}')
            
        # Generate interpolated points along the upper and lower surfaces
        upper_points, lower_points = [], []
        for t in parametric_points_up:
            # Place points on each edge using normalized parameter t
            upper_point  = upper_surface.positionAt(t).toTuple()
            upper_points.append(upper_point)
                
        for k in parametric_points_low:
            lower_point  = lower_surface.positionAt(k).toTuple()
            lower_points.append(lower_point)
            
        # Create coords
        upper_points = np.array(upper_points)
        lower_points = np.array(lower_points)
        
        # Insert Node ID  (From upper TE to lower TE)
        data_up  = np.insert(upper_points[:,:], 0, Node_ID_up, axis=1)    
        data_low = np.insert(lower_points[1:,:], 0, Node_ID_low,axis=1)
        
        data = np.append(data_up, data_low, axis=0)
        
        # Sort wrt ascending Node ID
        sorted_indices = np.argsort(data[:,0])
        data = data[sorted_indices]
        
        # Check for closing the TE[:,1:]*1e-03
        if close_TE is True:
            data = modify.close_TE_gap(data, Node_ID_one_surf, n=10)
        
        # Convert from mm to meter
        data[:,1:] = data[:,1:]*1e-03
        
        ### Append ALL data ###
        ALL_sections.append(data)
        # Cross sections append (Optional)
        ALL_cross_sections.append(cross_section.val())    
        
    ## Export all cross sections as a step file (Optional)
    all_sections_compound = cq.Compound.makeCompound(ALL_cross_sections)
    cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")        
    
    return ALL_sections
    
### Get Coordinates ###
ALL_sections = get_coords(stp_file, num_points, z_planes, remove_TE)

## Delete the second last span since it creates skew panels
if dist_spanwise == 'cosine_TIP':
    del ALL_sections[-2]

if close_TE is True:
    TE_property = '-closed_TE'
else:
    TE_property=''
    
# Write the all coordinates to ALL_sections txt file    
with open(f'{output_dir}/12x6_ClarkY_Blade-span_{dist_spanwise}-sec_{dist_airfoil}{TE_property}.pmt', 'w') as file:
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
        
with open(f'{output_dir}/Blade_points_check.txt', 'w') as file:
    file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
    for i in range(len(ALL_sections)):
        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')
        
###############################################################################################################################
end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds")