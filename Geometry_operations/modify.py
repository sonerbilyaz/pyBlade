import numpy as np

def remove_TE(edge_objects):
    # Get the midpoints of the edges
    count = 0
    for edge in edge_objects:
        count = count+1
        midpoint=np.array(edge.positionAt(0.5).toTuple())
        if count == 1:
            edges_midpoints = midpoint.reshape(1,3)
        if count != 1:
            edges_midpoints = np.vstack((edges_midpoints, midpoint))
            
    ## Select the TE whose midpoint has the highest x-coordinate ##
    TE_id = np.argmax(edges_midpoints[:,0])
    TE = edge_objects[TE_id]
    ## Get the TE coordinates
    p1=TE.positionAt(0).toTuple()                   ## Point coords at one vertex
    p2=TE.positionAt(1).toTuple()                   ## Point coords at the other vertex
    
    ## Amont TE coordinates, y-coordinate of the upper TE should be lower than the midpoint y-coordinate
    ## of the TE 
    if p2[1] < p1[1]:
        TE_upper = np.array(p2)
        TE_lower = np.array(p1)
        
    if p1[1] < p2[1]:
        TE_upper = np.array(p1)
        TE_lower = np.array(p2)
        
    del edge_objects[TE_id]
    filtered_edges = edge_objects
    
    return filtered_edges, TE_upper, TE_lower

def close_TE_gap(points, Node_IDs_upper_surface, n):
    
    ##################  Get the information from the section ##################
    num_points = points.shape[0]
    # First point in the sectional points is "TE_upper", and the last one is "TE_lower"
    TE_upper = points[0,1:]
    TE_lower = points[-1,1:]
    
    TE_midpoint = (TE_lower + TE_lower)/2
    TE_length = np.linalg.norm(TE_upper - TE_lower)
    
    # Last point in the upper surface is LE. Get this ID and find the coordinates at this ID in "points" array
    LE_ID = Node_IDs_upper_surface[0]
    LE_coordinates = points[points[:,0] == LE_ID][:,1:]

    # Get the chord vector and the chord length#
    chord_vector = LE_coordinates - TE_midpoint
    chord_length = np.linalg.norm(chord_vector)
    
    # Calculate the angle between the chord_vector and the -x axis
    x_axis = np.array([-1, 0, 0])
    dot_product = np.dot(chord_vector, x_axis)
    
    angle_radians = np.arccos(dot_product / chord_length)[0]
    
    ##########################  Modify the section ############################
    """" Move to origin """
    # Move points to the origin (LE is at the origin)
    LE_coordinates_shaped = LE_coordinates[0] * np.ones((num_points,1))
    points_origin = points[:,1:] - LE_coordinates_shaped
    
    """" Rotate """
    # 2D Rotation matrix for -z axis rotation N_wake=2-N_par=500-n_points=43-cosine_LE
    Rz_neg = np.array([
        [np.cos(angle_radians), np.sin(angle_radians)],
        [-np.sin(angle_radians),  np.cos(angle_radians)]
    ])
    # Rotate the points (only x and y coordinates)
    points_origin_rotated_xy = np.dot(points_origin[:,0:2], Rz_neg.T)
    
    """ Scale """
    # Scale the points to get the chord line between 0 and 1
    points_origin_rotated_xy_scaled = points_origin_rotated_xy / chord_length
    # Make sure that TE x coordinates are 1
    points_origin_rotated_xy_scaled[0,0], points_origin_rotated_xy_scaled[-1,0] = 1, 1
    
    """ Close TE """
    # Get the upper and lower surfaces coordinates
    upper_surface = points_origin_rotated_xy_scaled[0:LE_ID,:]    
    lower_surface = points_origin_rotated_xy_scaled[LE_ID:,:]
    
    ## Close the TE with the equation
    upper_surface[:,1] = upper_surface[:,1] + upper_surface[:,0]**n*((TE_length/chord_length)/2)
    lower_surface[:,1] = lower_surface[:,1] - lower_surface[:,0]**n*((TE_length/chord_length)/2)
    points_new_xy = np.vstack((upper_surface, lower_surface))
    
    #### REVERSE PROCESS ####
    """ Scale BACK """
    points_new_xy_scaled = points_new_xy * chord_length
    """ Rotate BACK """
    # 2D Rotation matrix for +z axis rotation 
    Rz = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians),  np.cos(angle_radians)]
    ])
    points_new_xy_scaled_rotated = np.dot(points_new_xy_scaled, Rz.T)
    
    # Add Node ID and z-coordinate
    points_new_scaled_rotated = np.hstack((np.linspace(1, num_points, num_points).reshape(num_points,1), points_new_xy_scaled_rotated)) # Insert Node ID
    points_new_scaled_rotated = np.hstack((points_new_scaled_rotated, np.zeros((num_points,1)))) # Insert z-coordinate
    
    """ Move to Original position """
    points_new_scaled_rotated[:,1:] = points_new_scaled_rotated[:,1:] + LE_coordinates_shaped
    
    points_new = points_new_scaled_rotated 
    
    return points_new