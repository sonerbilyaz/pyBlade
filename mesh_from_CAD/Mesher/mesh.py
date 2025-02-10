import meshio
import numpy as np
from Geometry_operations.Rotation import Rotate

def generate_mesh(ALL_sections, n_blades, close_TE):
    
    """ ELEMENT ORDER: Start from lower left and go clockwise (Top view looking to suction surface)
    
    UL    UR              LL (lower left), numpy array (3,)
    o --- o               UL (upper left), numpy array (3,)
    \     \               UR (upper right), numpy array (3,)
    o --- o               LR (lower right), numpy array (3,)
    LL    LR
    
    """
    
    num_points = len(ALL_sections[0])
    
    ## Reconstruct the coordinates array such that every span vertically concetaneted after the last one
    for i, section in enumerate(ALL_sections):
        if i == 0:
            ALL_sections_one_array = section
        else:
            ALL_sections_one_array = np.vstack((ALL_sections_one_array, section))
    
    
    ## Build the connectivity array, which contains the point IDs. shape = (cell_no ,4)
    for i in range(len(ALL_sections_one_array)-num_points):
        ID = i+1  # Which ID are we parsing?? 
        
        ## Connectivity will not change for the OPEN_TE case and the intermediate points in closed_TE case
        if close_TE is False or ID%num_points!=0:
            p1_index = i                 # LL
            p2_index = i+1               # UL
            p3_index = i+1+num_points       # UR
            p4_index = i+num_points         # LR
            
        ## Connectivity will be different for the last point in each section
        ## if we close the TE and remove the last point !!
        if close_TE is True and ID%num_points==0:
            p1_index = i                        # LL
            p2_index = i-(num_points-1)       # UL
            p3_index = i+1                      # UR
            p4_index = i+1+(num_points-1)     # LR
            
        cell_connect = np.array([[p1_index, p2_index, p3_index, p4_index]])
        
        cell_connect_DUST = np.array([[p1_index+1, p2_index+1, p3_index+1, p4_index+1]])
        
        if i == 0:
            connectivity = cell_connect
            connectivity_DUST = cell_connect_DUST
            
        if i != 0:
            connectivity = np.vstack((connectivity, cell_connect))
            connectivity_DUST = np.vstack((connectivity_DUST, cell_connect_DUST))
            
    cell_points = ALL_sections_one_array[:,1:]
    
    """ #########     CREATE THE COORDINATES OF OTHER BLADES     ########## """
    
    ## Construct the points of the other blades by rotating the "cell_points" ##
    if n_blades != 1:
        # Construct the list which contains coordinates of each blade
        cell_points_ALL = [cell_points for i in range(n_blades)]
        
        # Rotation Angle #
        angle = 360/n_blades
        
        ## Rotate each element to get the coordinates of each blade ##
        for i in range(n_blades-1):
            # Rotate the points #
            rotated_points = Rotate(cell_points_ALL[i], angle, axis='y')
            # Update the new blade points
            cell_points_ALL[i+1] = rotated_points
    
    elif n_blades == 1:
        cell_points_ALL = [cell_points]
        
    ## Convert the cell_points_ALL array and connectivity array into a single array ##
    for i in range(len(cell_points_ALL)):
        if i == 0:
            combined_points = cell_points_ALL[i]
            combined_connectivity = connectivity
            combined_connectivity_DUST = connectivity_DUST
        if i != 0:
            combined_points = np.vstack((combined_points, cell_points_ALL[i]))
            combined_connectivity = np.vstack((combined_connectivity, connectivity + i*cell_points.shape[0]*np.ones(connectivity.shape)))
            combined_connectivity_DUST = np.vstack((combined_connectivity_DUST, connectivity_DUST + i*cell_points.shape[0]*np.ones(connectivity_DUST.shape)))
    
    """ ##############         CREATE MESH FOR EACH BLADE      #############"""
    # Generate the DUST mesh #
    combined_points_DUST = Rotate(combined_points, -90, axis='x')
    combined_points_DUST = Rotate(combined_points_DUST, -90, axis='z')
    
    cells = [('quad', combined_connectivity)]
    mesh_DUST= meshio.Mesh(combined_points_DUST, cells)
    
    # Generate the overall mesh #
    cells = [('quad', combined_connectivity)]
    mesh= meshio.Mesh(combined_points, cells)
    
    return mesh, mesh_DUST, combined_points_DUST, combined_connectivity_DUST
