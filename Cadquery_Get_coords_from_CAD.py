from IPython import get_ipython
get_ipython().run_line_magic('clear', '')

from timeit import default_timer as timer
start_time = timer()  # Start the timer

import cadquery as cq
import numpy as np
import os 

"""################################ INPUTS #######################################"""
# File paths and the STEP file
working_dir = 'Test_Cases/eVTOLUTION_VX4_Front_Propeller'
stp_file = 'Front_Propeller_Blade.stp'

output_dir = f'{working_dir}/output'
############################ Geometry Parameters ##############################
surf_type = 'propeller'     # Surface type (wing/propeller)
number_of_blades = 5        # Number of blades
remove_TE = True            # Should we remove TE ??

rotation_center = [0, 0, 0]
rotation_axis = [0, -1, 0]

############################ Panel Parameters ##################################
## Airfoil Section ##
# Number of points to generate (upper and lower surf separately)
num_points = 21
dist_airfoil = 'cosine_LE'

## Spanwise Cutting Planes ##
spanwise_panel_num= 32
z_min, z_max = 26.2, 149 ## in mm
dist_spanwise = 'cosine_TIP'

"""###############################################################################"""
# Create output directory
if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)
    
# Normalizing function when we shift clustering location of the points between (0 and 1)
def normalize(data, min_val, max_val):
    return min_val + (data - data.min()) * (max_val - min_val) / (data.max() - data.min())

def points_distribution(num_points, dist = dist_airfoil): 
    if dist == 'cosine_LE-TE':
        # Generate evenly spaced angles between 0 and 90 deg
        angles  = np.linspace(np.radians(0), np.radians(90), num_points//2)
        # Divide the first half (between 0-0.5) into points which cluster towards both ends (0 and 0.5)
        x_half_1, x_half_2= (1-0.5*np.sin(angles)), (0.5*np.sin(angles))

        # Remove repetitive point (midpoint = 0.5)
        x_half_1 = np.delete(np.sort(x_half_1),-1)
        
        # Connect the halves
        points = np.sort(np.hstack((x_half_1,x_half_2)))
        
        # if num_points is NOT an even number
        if num_points/2 != num_points//2:
            Node_ID_one_surf = np.linspace(num_points-2, 1, num_points-2)
            
        # if num_points is an even  number
        if num_points/2 == num_points//2:
            Node_ID_one_surf = np.linspace(num_points-1, 1, num_points-1)
            
    if dist == 'cosine_LE':    
        # Generate evenly spaced angles between 0 and 90 deg
        angles  = np.linspace(np.radians(0), np.radians(90), num_points)
        # Get the distribution
        points = 1-np.sin(angles)
        # Sort in ascending order
        points = np.sort(points)
        
        Node_ID_one_surf = np.linspace(num_points, 1, num_points)
    if dist == 'linear':
        # Remain linear
        points = np.linspace(1, 0, num_points)
        
        Node_ID_one_surf = np.linspace(num_points, 1, num_points)
    
    return points, Node_ID_one_surf

def identify_edges(filtered_edges):
    if len(filtered_edges) != 2:
        raise ValueError("There should be exactly two edges to merge.")
    
    # Find shared vertex
    vertices = [cq.Vector(v.toTuple()) for edge in filtered_edges for v in edge.vertices()]
    common_vertex = next((v for v in vertices if vertices.count(v) > 1), None)
    
    if not common_vertex:
        raise ValueError("No common vertex found between the two edges.")
    
    # Identify LE coordinates
    x,y,z = common_vertex.toTuple()
    LE_coordinates = np.array([x,y,z])
    
    # Identify non-common vertices (TE) for each edge
    edges_with_non_common_vertices = []
    for edge in filtered_edges:
        for v in edge.vertices():
            if cq.Vector(v.toTuple()) != common_vertex:
                edges_with_non_common_vertices.append([edge, cq.Vector(v.toTuple())])
                
    # Compare the y-coordinates of the non-common vertices (TE). Lower-y coord is UPPER surface
    if edges_with_non_common_vertices[0][1].y < edges_with_non_common_vertices[1][1].y:
        upper_surface, TE_upper_coordinates = edges_with_non_common_vertices[0][0], np.array([edges_with_non_common_vertices[0][1].toTuple()])
        lower_surface, TE_lower_coordinates = edges_with_non_common_vertices[1][0], np.array([edges_with_non_common_vertices[1][1].toTuple()])
    elif edges_with_non_common_vertices[0][1].y > edges_with_non_common_vertices[1][1].y:
        upper_surface, TE_upper_coordinates = edges_with_non_common_vertices[1][0], np.array([edges_with_non_common_vertices[1][1].toTuple()])
        lower_surface, TE_lower_coordinates = edges_with_non_common_vertices[0][0], np.array([edges_with_non_common_vertices[0][1].toTuple()])
    
    return upper_surface, lower_surface, LE_coordinates, TE_upper_coordinates, TE_lower_coordinates

def spanwise_disribution(z_min, z_max, spanwise_panel_num, r_R = 0.8, dist = dist_spanwise):
    if dist == "cosine":
        # Generate evenly spaced angles between 0 and π
        angles = np.linspace(0, np.pi, spanwise_panel_num)
        # Apply cosine function and normalize to [0, 1]
        points = 0.5 * (1 - np.cos(angles))
        # Map points to [z_min, z_max]
        points = z_min + points * (z_max - z_min)
        return points
    
    if dist == "cosine_TIP":
        angles  = np.linspace(np.radians(0), np.radians(90), 2*int((1-r_R) * spanwise_panel_num))
        
        # Uniform part towards root. Remove extra end point
        x1 = np.linspace(0, r_R, int(r_R * spanwise_panel_num))
        x1 = np.delete(x1, -1)
        # cosine TIP part
        x2 = normalize(np.sin(angles), r_R, 1)
        # Connect them
        points = np.hstack((x1,x2))
        # Map points to [z_min, z_max]
        points = z_min + points * (z_max - z_min)
        return points
    
    elif dist == "uniform":
        return np.linspace(z_min, z_max, spanwise_panel_num)      # Linear Distribution

z_planes = spanwise_disribution(z_min, z_max, spanwise_panel_num)

def remove_TE_from_CAD(cross_sec_edge_objects):
    # Calculate lengths of edges
    edge_lengths = [(edge, edge.Length()) for edge in cross_sec_edge_objects]
    # Find the shortest edge (TE)
    shortest_edge = min(edge_lengths, key=lambda x: x[1])[0]
    
    # Filter out the shortest edge
    filtered_edges = [edge for edge, length in edge_lengths if edge != shortest_edge]
    
    return filtered_edges
    
         
def get_coords(stp_file, num_points, z_planes, remove_TE):    
    # Load the .stp file
    blade = cq.importers.importStep(f'{working_dir}/{stp_file}')
    ALL_cross_sections=[]
    ALL_sections=[]
    
    
    for z in z_planes:
        # Create a section of the blade using the defined z-coord
        cross_section=blade.section(height=z)
        # Append the ALL cross sections and export (Optional)
        ALL_cross_sections.append(cross_section.val())        
        # Get all edges and calculate their lengths
        edges = cross_section.edges()
        # Access the edge objects from the Workplane
        edge_objects = edges.objects
        
        # Should we remove TE ??
        if remove_TE is True:
            filtered_edges = remove_TE_from_CAD(edge_objects)
            # Check that there are exactly 2 edges remaining
            if len(filtered_edges) != 2:
                print("ERROR: Cross section does not have exactly 2 edges after filtering.")
                break
            
        # If we don't remove, check the detected edge numbers
        if remove_TE is False and len(edge_objects) != 2:
            print(f'ERROR: There are not 2 edges detected at cross section z={z}mm !!')
            break
        if remove_TE is False and len(edge_objects) == 2:
            filtered_edges = edge_objects
            
        # Identify the edges
        upper_surface, lower_surface, LE_coords, TE_upper_coords, TE_lower_coords = identify_edges(filtered_edges)
        
        ## Generate parametric points between 0 and 1
        parametric_points_up, Node_ID_one_surf = points_distribution(num_points)
        parametric_points_low, Node_ID_one_surf = points_distribution(num_points)
        
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
        ALL_sections.append(data)
        
        
        # ALL_upper_suface.append(upper_surface)
        # ALL_parametric.append([parametric_points_up, parametric_points_low])
        # ALL_LE_coords.append(LE_coords)
        # ALL_TE_upper_coords.append(TE_upper_coords)
        # ALL_TE_lower_coords.append(TE_lower_coords)
        # ALL_upper_points.append(upper_points)
        
    ## Export all cross sections as a step file (Optional)
    all_sections_compound = cq.Compound.makeCompound(ALL_cross_sections)
    cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")        
    
    # ## Convert ALL coordinates from mm to meter
    # for i in range(len(ALL_sections)):
    #     ALL_sections[i][:,1:] = ALL_sections[i][:,1:]*1e-03
    
    return ALL_sections, len(Node_ID_one_surf)
    
### Get Coordinates ###
ALL_sections, Nodes_one_surf = get_coords(stp_file, num_points, z_planes, remove_TE)

# Write the all coordinates to ALL_sections txt file    
with open(f'{output_dir}/VX4_Front_1_Blade-span_{dist_spanwise}-sec_{dist_airfoil}.pmt', 'w') as file:
    file.write('######## Panel parameters ########\n')
    file.write('type=' + surf_type + '\n')
    file.write('n_blades=' + str(number_of_blades) + '\n\n')
    file.write('rotation_center=' + str(rotation_center) + '\n')
    file.write('rotation_axis=' + str(rotation_axis) + '\n')
    
    file.write('n_span_all=' + str(len(z_planes)) + '\n')
    file.write('n_points=' + str(2*Nodes_one_surf-1) + '\n')
    file.write('######## End of parameters ########\n')
    # file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
    for i in range(len(ALL_sections)):
        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.5f','%.5f','%.5f'], comments='')
        
with open(f'{output_dir}/Blade_points_check.txt', 'w') as file:
    file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
    for i in range(len(ALL_sections)):
        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.5f','%.5f','%.5f'], comments='')

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