import meshio
import numpy as np

def generate_mesh(ALL_sections):
    span = len(ALL_sections)
    num_points_sec = len(ALL_sections[0])
    
    """ ELEMENT ORDER: Start from lower left and go clockwise 
    
    UL    UR              LL (lower left), numpy array (3,)
    o --- o               UL (upper left), numpy array (3,)
    \     \               UR (upper right), numpy array (3,)
    o --- o               LR (lower right), numpy array (3,)
    LL    LR
    
    """
    
    ## Rearrange the point coordinates in a suitable format ##
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
    
    # Convert to an integer array
    order = order.astype(int)
    ### CREATE MESH ###
    cells = [('quad',order)]
    mesh = meshio.Mesh(cell_points, cells)
    
    return mesh