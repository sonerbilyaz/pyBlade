import meshio
import numpy as np

def generate_mesh(ALL_sections, n_blades):
    span = len(ALL_sections)
    num_points_sec = len(ALL_sections[0])
    
    """ ELEMENT ORDER: Start from lower left and go clockwise 
    
    UL    UR              LL (lower left), numpy array (3,)
    o --- o               UL (upper left), numpy array (3,)
    \     \               UR (upper right), numpy array (3,)
    o --- o               LR (lower right), numpy array (3,)
    LL    LR
    
    """
    
    ## Rearrange the point coordinates of a 1 BLADE in a suitable format ##
    count = 0
    for ID in range(num_points_sec-1):
        for i in range(span-1):
            count = count +1
            
            ### LL-UL-UR-LR ####
            p1 = ALL_sections[i][ID,1:]         # LL
            p2 = ALL_sections[i][ID+1,1:]       # UL    
            p3 = ALL_sections[i+1][ID+1,1:]     # UR    
            p4 = ALL_sections[i+1][ID,1:]       # LR    
            
            points = np.array([p1,p2,p3,p4])
            
            if count == 1:
                cell_points = points
                
            if count != 1:
                cell_points = np.vstack((cell_points, points))
                
    ## Construct the array which contains the indices of these ordered coordinates ##
    indices = np.linspace(0, len(cell_points)-1, len(cell_points))
    
    # Rearrange #
    for i in range(0, len(cell_points)-3, 4):
        if i == 0:
            order = indices[i:i+4].reshape(1,4)
        if i != 0:
            order = np.vstack((order, indices[i:i+4].reshape(1,4)))
    
    # Convert to an integer array. Each row in this array contains the indices of the points which are ordered to construct one cell
    order = order.astype(int)
    
    """   #####              CREATE OTHER BLADES        #######             """
    
    ## Construct the points of the other blades by rotating the "cell_points" ##

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
    
    """ ##############         CREATE MESH FOR EACH BLADE      #############"""
    # Combine the cell_points_ALL and order into a single array #
    for i in range(len(cell_points_ALL)):
        if i == 0:
            combined_points = cell_points_ALL[i]
            combined_order = order
        if i != 0:
            combined_points = np.vstack((combined_points, cell_points_ALL[i]))
            combined_order = np.vstack((combined_order, order + i*cell_points.shape[0]*np.ones(order.shape)))
        
    # Generate the mesh for the propeller
    cells = [('quad', combined_order)]
    propeller = meshio.Mesh(combined_points, cells)
    
    return propeller