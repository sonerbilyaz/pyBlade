import numpy as np
from . import extract
from .Rotation import Rotate

def remove_TE(edge_objects):
    """
    Removes the TE edge object from all of the edges
        
        Parameters
        ----------
        edge_objects : List
            All edge objects of a section (includes TE)

        Returns
        -------
        filtered_edges: List
            All edge objects of a section (WITHOUT TE!!)
    """

    ## Get the TE
    TE, _, _ = extract.get_TE(edge_objects)

    ## Remove the TE from the edges
    for i, edge in enumerate(edge_objects):
        if edge == TE:
            del edge_objects[i]

    filtered_edges = edge_objects
    
    return filtered_edges
    
def pitch_increase(points, pitch_increment):
    points = np.concatenate([points, np.zeros((points.shape[0],1))], axis=1)
    points_new = Rotate(points, pitch_increment, 'z')

    return points_new[:,:2]

def sweep_increase(points, sweep_deg, z_span, r_root):
    
    sweep = np.tan(np.radians(sweep_deg)) * (z_span - r_root)

    points[:,1] = points[:,1] + sweep

    return points


def modify_airfoil(points, z_span, Node_IDs_upper_surface, coll_pitch, scale, sweep, r_root, n=10, x_0 = 0.7):
    
    # -----------------------  Get the information from the section -----------------------
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
    
    quart_chord_coords = LE_coordinates - chord_vector*0.25

    # Calculate the angle between the chord_vector and the -x axis
    x_axis = np.array([-1, 0, 0])
    dot_product = np.dot(chord_vector, x_axis)
    
    angle_radians = np.arccos(dot_product / chord_length)[0]
    # ---------------------------------------------------------------------------------------

    """ ##################    MOVING TO THE ORIGIN (NORMALIZE)    ################## """
    # Move points to the origin (0.25*chord is at the origin)
    points_origin = points[:,1:] - quart_chord_coords * np.ones_like(points[:,1:])
    
    ### Rotate ###
    # 2D Rotation matrix for -z axis rotation
    Rz_neg = np.array([
        [np.cos(angle_radians), np.sin(angle_radians)],
        [-np.sin(angle_radians),  np.cos(angle_radians)]
    ])
    # Rotate the points (only x and y coordinates)
    points_origin_rotated_xy = np.dot(points_origin[:,:2], Rz_neg.T)

    ### Scale ###
    # Scale the points to get the chord line between 0 and 1
    points_origin_xy = points_origin_rotated_xy / chord_length
    # Make sure that TE x coordinates are 0.75 (origin is 0.25*chord!)
    points_origin_xy[0,0], points_origin_xy[-1,0] = 0.75, 0.75
    """ ############################################################################  """

    #######     1) CLOSE TE ??      #######
    # Get the upper and lower surfaces coordinates
    upper_surface = points_origin_xy[0:LE_ID,:]    
    lower_surface = points_origin_xy[LE_ID:,:]
    
    ## Close the TE with the equation
    ## TE closure function (Origin should be at the LE !!!!)
    x_up, x_low = upper_surface[:,0] + 0.25, lower_surface[:,0] + 0.25
    w_x_up = np.where(x_up > x_0, x_up**n, 0)
    w_x_low = np.where(x_low > x_0, x_low**n, 0)

    upper_surface[:,1] = upper_surface[:,1] + w_x_up*(TE_length/chord_length)/2
    lower_surface[:,1] = lower_surface[:,1] - w_x_low*(TE_length/chord_length)/2
    points_origin_xy = np.vstack((upper_surface, lower_surface))
    
    #######     2) COLL PITCH        #######
    points_origin_xy = pitch_increase(points_origin_xy, coll_pitch)

    #######     3) CHORD SCALE       #######
    points_origin_xy = points_origin_xy * scale

    
    """ #########################  REVERSE PROCESS #########################    """ 
    ### Scale BACK ###
    points_new_xy = points_origin_xy * chord_length
    
    ### Rotate BACK ###
    # 2D Rotation matrix for +z axis rotation 
    Rz = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians),  np.cos(angle_radians)]
    ])
    
    points_new_xy = np.dot(points_new_xy, Rz.T)

    # Add Node ID and z-coordinate
    points_new = np.hstack((np.linspace(1, num_points, num_points).reshape(num_points,1), points_new_xy)) # Insert Node ID
    points_new = np.hstack((points_new, np.zeros((num_points,1)))) # Insert z-coordinate
    
    ### Move to Original position ###
    points_new[:,1:] = points_new[:,1:] + quart_chord_coords * np.ones_like(points[:,1:])

    #######     4) SWEEP             #######
    points_new = sweep_increase(points_new, sweep, z_span, r_root)

    return points_new, np.rad2deg(angle_radians), chord_length

def change_span(points, Node_IDs_upper_surface, twist_local, taper_local):

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
    ### Move to origin ###
    # Move points to the origin (LE is at the origin)
    LE_coordinates_shaped = LE_coordinates[0] * np.ones((num_points,1))
    points_origin = points[:,1:] - LE_coordinates_shaped
    
    ### Rotate ###
    # 2D Rotation matrix for -z axis rotation
    Rz_neg = np.array([
        [np.cos(angle_radians), np.sin(angle_radians)],
        [-np.sin(angle_radians),  np.cos(angle_radians)]
    ])
    # Rotate the points (only x and y coordinates)
    points_origin_rotated_xy = np.dot(points_origin[:,0:2], Rz_neg.T)
    
    # Center the section at the QUARTER chord location
    ex = np.array([[1,0]])       # Unit vector in x-dir
    points_origin_rotated_xy = points_origin_rotated_xy - ex*chord_length*0.25
    
    ### Scale  ###
    points_new_xy_scaled = points_origin_rotated_xy * taper_local

    ### Add twist ###
    # 2D Rotation matrix for +z axis rotation 
    Rz = np.array([
        [np.cos(np.radians(twist_local)), -np.sin(np.radians(twist_local))],
        [np.sin(np.radians(twist_local)),  np.cos(np.radians(twist_local))]
    ])

    points_new_xy_scaled_twisted = np.dot(points_new_xy_scaled, Rz.T)

    """ #########################  REVERSE PROCESS #########################    """ 
    ### Move LE to Origin ###
    points_new_xy = points_new_xy_scaled_twisted + ex*chord_length*0.25
    
    ### Rotate BACK ###
    # 2D Rotation matrix for +z axis rotation 
    Rz = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians),  np.cos(angle_radians)]
    ])
    
    points_new_origin = np.dot(points_new_xy, Rz.T)

    # Add Node ID and z-coordinate
    points_new_scaled_rotated = np.hstack((np.linspace(1, num_points, num_points).reshape(num_points,1), points_new_origin)) # Insert Node ID
    points_new_scaled_rotated = np.hstack((points_new_scaled_rotated, np.zeros((num_points,1)))) # Insert z-coordinate
    
    ### Move to Original position ###
    points_new_scaled_rotated[:,1:] = points_new_scaled_rotated[:,1:] + LE_coordinates_shaped
    
    ## Make sure that the upper TE node and lower TE node are the same (1st and last point will be the same, coincides on top of each other)
    points_new = points_new_scaled_rotated.copy() 
    return points_new
