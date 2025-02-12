import meshio
import numpy as np
from Geometry_operations.Rotation import Rotate

def generate_mesh(ALL_sections, ALL_sections_DUST,  n_blades, close_TE):
    
    """ ELEMENT ORDER: Start from lower left and go clockwise (Top view looking to suction surface)
    
    UL    UR              LL (lower left), numpy array (3,)
    o --- o               UL (upper left), numpy array (3,)
    \     \               UR (upper right), numpy array (3,)
    o --- o               LR (lower right), numpy array (3,)
    LL    LR
    
    """
    
    num_points = len(ALL_sections[0])
    num_points_DUST = len(ALL_sections_DUST[0])
    
    ## Reconstruct the coordinates array such that every span vertically concetaneted after the last one
    count = 0
    for section, section_DUST in zip(ALL_sections, ALL_sections_DUST):
        count = count + 1
        
        if count == 1:
            ALL_sections_one_array = section
            ALL_sections_one_array_DUST = section_DUST
        else:
            ALL_sections_one_array = np.vstack((ALL_sections_one_array, section))
            ALL_sections_one_array_DUST = np.vstack((ALL_sections_one_array_DUST, section_DUST))
        
    #####   Build the connectivity array contains the point IDs. shape = (cell_no ,4)     ######
    for i,j in zip(range(len(ALL_sections_one_array)-num_points),range(len(ALL_sections_one_array_DUST)-num_points_DUST)):        
        ###### ---- NVLM Build the connectivity array ---- #######
        p1_index = i                                # LL
        p2_index = i+1                              # UL
        p3_index = i+1+num_points                   # UR
        p4_index = i+num_points                     # LR
        
        ###### ---- DUST Build the connectivity array ---- #######
        ID = j+1
        # Connectivity will not change for the intermediate points
        if close_TE is False or ID%num_points_DUST!=0:
            p1_index_DUST = j                       # LL
            p2_index_DUST = j+1                     # UL
            p3_index_DUST = j+1+num_points_DUST     # UR
            p4_index_DUST = j+num_points_DUST       # LR
        
        # Connectivity will be different for the last TE point !!!!!! 
        if close_TE is True and ID%num_points_DUST==0:
            p1_index_DUST = j                      # LL
            p2_index_DUST = j-(num_points_DUST-1)       # UL
            p3_index_DUST = j+1                    # UR
            p4_index_DUST = j+1+(num_points_DUST-1)     # LR

        
        cell_connect = np.array([[p4_index, p3_index, p2_index, p1_index]])
        cell_connect_DUST = np.array([[p4_index_DUST, p3_index_DUST, p2_index_DUST, p1_index_DUST]])
        
        if i == 0:
            connectivity = cell_connect
        if i != 0:
            connectivity = np.vstack((connectivity, cell_connect))
            
        if j == 0:
            connectivity_DUST = cell_connect_DUST
        if j != 0:
            connectivity_DUST = np.vstack((connectivity_DUST, cell_connect_DUST))
            
    cell_points = ALL_sections_one_array[:,1:]
    cell_points_DUST = ALL_sections_one_array_DUST[:,1:]
    
    """ ##############         CREATE MESH       #############"""
    # Generate the DUST mesh #
    cell_points_DUST = Rotate(cell_points_DUST, -90, axis='x')
    cell_points_DUST = Rotate(cell_points_DUST, -90, axis='z')
    
    cells = [('quad', connectivity_DUST)]
    mesh_DUST= meshio.Mesh(cell_points_DUST, cells)
    
    # Generate the overall mesh #
    cells = [('quad', connectivity)]
    mesh= meshio.Mesh(cell_points, cells)
    
    return mesh, connectivity, mesh_DUST, connectivity_DUST, cell_points_DUST
