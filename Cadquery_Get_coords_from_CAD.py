from IPython import get_ipython
get_ipython().run_line_magic('clear', '')

from timeit import default_timer as timer
start_time = timer()  # Start the timer

from Distributions.distributions import airfoil_distribution as airfoil_points
from Distributions.distributions import spanwise_disribution as spanwise_planes

import Geometry_operations.modify as modify
import Geometry_operations.extract_info as extract

import cadquery as cq
import numpy as np, os 

"""################################ INPUTS #######################################"""
# File paths and the STEP file
working_dir = 'Test_Cases/eVTOLUTION_VX4_Front_Propeller'
stp_file = f'{working_dir}/Front_Propeller_Blade.stp'

############################ Geometry Parameters ##############################
surf_type = 'propeller'     # Surface type (wing/propeller)
number_of_blades = 5        # Number of blades
remove_TE = True            # Should we remove TE ??
close_TE = True             # Should we close the TE gap ??

rotation_center = [0, 0, 0]
rotation_axis = [0, -1, 0]

############################ Panel Parameters ##################################
""" Airfoil Section """
# Number of points to generate (upper and lower surf separately)
num_points = 20
dist_airfoil = 'cosine_LE'

""" Spanwise Cutting Planes """
spanwise_panel_num= 32
z_min, z_max = 26.2, 149            ## in mm
dist_spanwise = 'cosine_TIP'
r_R = 0.8                           ## Start span location of the cosine_TIP 
  
"""###############################################################################"""


output_dir = f'{working_dir}/output'
# Create output directory if there is no #
if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

### Create parametric points for airfoil distribution ###
parametric_points, Node_ID_one_surf = airfoil_points(num_points, dist_airfoil)

### Create spanwise planes ###
z_planes = spanwise_planes(z_min, z_max, spanwise_panel_num, dist_spanwise, r_R)

## Get coordinates ##
def get_coords(stp_file, num_points, z_planes, remove_TE):    
    # Load the .stp file
    blade = cq.importers.importStep(stp_file)
    
    ALL_sections = []
    ALL_sections_closed_TE = []
    ALL_cross_sections = []
    
    for z in z_planes:
        # Create a section of the blade using the defined z-coord
        cross_section=blade.section(height=z)
        # Get all edges from the cross section
        edges = cross_section.edges()
        # Access the edge objects from the Workplane
        edge_objects = edges.objects
        
        # Should we remove TE ??
        if remove_TE is True:
            filtered_edges = modify.remove_TE_from_CAD(edge_objects)
            # Check that there are exactly 2 edges remaining
            if len(filtered_edges) != 2:
                print("ERROR: Cross section does not have exactly 2 edges after filtering.")
                break
            
        # If we don't remove, check the detected edge numbers
        if remove_TE is False and len(edge_objects) != 3:
            print(f'ERROR: There are not 3 edges (upper surface + lower surface + TE) detected at cross section z={z}mm !!')
            break
        if remove_TE is False and len(edge_objects) == 3:
            filtered_edges = edge_objects
            
        # Identify the edges
        upper_surface, lower_surface, LE_coords = extract.identify_edges(filtered_edges)
        
        ## Assign the parametric points between 0 and 1 to the upper and lower surfaces
        parametric_points_up = parametric_points 
        parametric_points_low = parametric_points
        
        # Check it starts generating points from the starting vertex (LE)
        if (np.round(upper_surface.positionAt(0).toTuple(),decimals=2) != np.round(LE_coords,decimals=2)).any() and dist_airfoil == 'cosine_LE':
            # If the orientation is wrong, reverse it to correct ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
            parametric_points_up = 1 - parametric_points_up
            
        if (np.round(lower_surface.positionAt(0).toTuple(),decimals=2) != np.round(LE_coords,decimals=2)).any() and dist_airfoil == 'cosine_LE':
            # If the orientation is wrong, reverse it to correct ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
            parametric_points_low = 1- parametric_points_low
        
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
        data_up  = np.insert(upper_points[:,:], 0, Node_ID_one_surf, axis=1)    
        data_low = np.insert(lower_points[1:,:], 0, (Node_ID_one_surf +len(Node_ID_one_surf)-1)[:-1][::-1],axis=1)
        
        data = np.append(data_up, data_low, axis=0)
        
        # Sort wrt ascending Node ID
        sorted_indices = np.argsort(data[:,0])
        data = data[sorted_indices]
        
        ## Closed TE ##
        data_closed_TE = modify.close_TE_gap(data, Node_ID_one_surf, n=10)
        
        ### Append ALL data ###
        ALL_sections.append(data)
        ALL_sections_closed_TE.append(data_closed_TE)
        # Cross sections export (Optional)
        ALL_cross_sections.append(cross_section.val())    
        
        
    # ## Export all cross sections as a step file (Optional)
    # all_sections_compound = cq.Compound.makeCompound(ALL_cross_sections)
    # cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")        
    
    # ## Convert ALL coordinates from mm to meter
    # for i in range(len(ALL_sections)):
    #     ALL_sections[i][:,1:] = ALL_sections[i][:,1:]*1e-03
    
    return ALL_sections, ALL_sections_closed_TE
    
### Get Coordinates ###
ALL_sections, ALL_sections_closed_TE = get_coords(stp_file, num_points, z_planes, remove_TE)

# Write the all coordinates to ALL_sections txt file    
with open(f'{output_dir}/VX4_Front_1_Blade-span_{dist_spanwise}-sec_{dist_airfoil}.pmt', 'w') as file:
    file.write('######## Panel parameters ########\n')
    file.write('type=' + surf_type + '\n')
    file.write('n_blades=' + str(number_of_blades) + '\n\n')
    file.write('rotation_center=' + str(rotation_center) + '\n')
    file.write('rotation_axis=' + str(rotation_axis) + '\n')
    
    file.write('n_span_all=' + str(len(z_planes)) + '\n')
    file.write('n_points=' + str(len(ALL_sections[0][:,0])) + '\n')
    file.write('######## End of parameters ########\n')
    
    for i in range(len(ALL_sections)):
        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.5f','%.5f','%.5f'], comments='')
        
with open(f'{output_dir}/Blade_points_check.txt', 'w') as file:
    file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
    for i in range(len(ALL_sections)):
        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.5f','%.5f','%.5f'], comments='')
        
with open(f'{output_dir}/Blade_points_check_closed_TE.txt', 'w') as file:
    file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
    for i in range(len(ALL_sections_closed_TE)):
        np.savetxt(file, ALL_sections_closed_TE[i], delimiter='\t', fmt=['%.0f','%.5f','%.5f','%.5f'], comments='')

# A_points=np.array([[]])
# for t in ALL_parametric[9][0]:
#     p = ALL_upper_suface[9].positionAt(t).toTuple()    
#     A_points = np.append(A_points, np.array([p]))

# A_data=np.insert(ALL_upper_points[9][1:,:], 0, np.linspace(num_points-1, 1, num_points-1),axis=1)    


# # Sort wrt ascending Node ID
# sorted_indices = np.argsort(A_data[:,0])
# A_data_sorted = A_data[sorted_indices]

###############################################################################################################################
end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds")