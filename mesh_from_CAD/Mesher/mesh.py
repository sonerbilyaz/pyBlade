import meshio
import numpy as np

def generate_mesh(ALL_sections, n_blades, close_TE):
    
    """ ELEMENT ORDER: Start from lower left and go clockwise (Top view looking to suction surface)
    
    UL    UR              LL (lower left), numpy array (3,)
    o --- o               UL (upper left), numpy array (3,)
    \     \               UR (upper right), numpy array (3,)
    o --- o               LR (lower right), numpy array (3,)
    LL    LR
    
    """
    
    """ Next, build connectivity !! After that, make sure that 
    "close_TE_gap" function removes the extra point after closing the TE !!"""
    
    
    ## Reconstruct the coordinates array such that every span vertically concetaneted after the last one
    ## (ONE POINT COORDINATES SHOULD NOT BE REPEATED AGAIN!)
    
    for i, section in enumerate(ALL_sections):
        if i == 0:
            ALL_sections_one_array = section
            ID_last = int(section[-1,0])
        else:
            section[:,0] = ALL_sections[i-1][:,0] + ID_last*np.ones((ID_last)) 
            ALL_sections_one_array = np.vstack((ALL_sections_one_array, section))
    
    
    ## Build the connectivity array, which contains the point IDs. shape = (cell_no ,4)
    for i in range(len(ALL_sections_one_array)-len(section)-1):
        ID = i+1  # Which ID are we parsing?? (1st column of ALL_sections_one_array)
        
        ## Connectivity will not change for the OPEN_TE case and the intermediate points in closed_TE case
        if close_TE is False or ID%len(section)!=0:
            p1_index = i                 # LL
            p2_index = i+1               # UL
            p3_index = i+1+ID_last       # UR
            p4_index = i+ID_last         # LR
            
        ## Connectivity will be different for the last point in each section
        ## if we close the TE and remove the last point !!
        if close_TE is True and ID%len(section)==0:
            p1_index = i                        # LL
            p2_index = i-(len(section)-1)       # UL
            p3_index = i+1                      # UR
            p4_index = i+1+(len(section)-1)     # LR
            
        cell_connect = np.array([[p1_index, p2_index, p3_index, p4_index]])
        
        ## The connectivities for DUST will be +1 of Python indices !!!  
        ## (They will be point IDs, 1st column of ALL_sections_one_array)
        
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
        
        # Rotation Matrix #
        angle = np.radians(360/n_blades)
        
        Ry = np.array([
        [np.cos(angle), 0, np.sin(angle)],
        [0,             1, 0            ],
        [-np.sin(angle), 0, np.cos(angle)]
        ])
        
        ## Rotate each element to get the coordinates of each blade ##
        for i in range(n_blades-1):
            # Rotate the points #
            rotated_points = np.dot(cell_points_ALL[i], Ry.T)
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
    # Reverse the connectivity order (CCW) (If necessary)
    combined_connectivity = combined_connectivity[:,::-1]
    combined_connectivity_DUST = combined_connectivity_DUST[:,::-1]
    
    # Generate the overall mesh #
    cells = [('quad', combined_connectivity)]
    mesh= meshio.Mesh(combined_points, cells)
    
    return mesh, combined_points, combined_connectivity_DUST
