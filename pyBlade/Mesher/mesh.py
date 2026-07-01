import meshio
import numpy as np
from ..Geometry_operations.Rotation import Rotate

def generate_mesh(ALL_sections, ALL_sections_DUST, close_TE):
    
    """ CCW ORDER: Start from lower left and go clockwise 
    (Looking from + to - surface normal !!)
    
    1     4               1 (upper left), numpy array (3,)
    o --- o               2 (lower left), numpy array (3,)
    |     |               3 (lower right), numpy array (3,)
    o --- o               4 (upper right), numpy array (3,)
    2     3
    
    """
    n_span = len(ALL_sections)
    n_chord = ALL_sections[0].shape[0]
    n_chord_DUST = ALL_sections_DUST[0].shape[0]

    blade = np.concatenate([section for section in ALL_sections], axis = 0)
    blade_DUST = np.concatenate([section for section in ALL_sections_DUST], axis = 0)

    ## (n_span, n_chord, 4) ==> 1st col Node ID
    blade_grid = blade.reshape((n_span, n_chord, 4))
    blade_DUST_grid = blade_DUST.reshape((n_span, n_chord_DUST, 4))

    connectivity = []
    connectivity_DUST = []

    for k in range(n_span-1):
        for i, i_DUST in zip(range(n_chord-1), range(n_chord_DUST)):
            
            ### ----------------------------    NVLM    ---------------------------- #####
            ## These indices will be the row rumbers of the point in blade array ==> (N,4) 
            ## Return the row number of the element in "blade" where the required one in "blade_grid" equals
            p1_index = np.flatnonzero((blade == blade_grid[k,i]).all(1))[0]         ## UL (1)   
            p2_index = np.flatnonzero((blade == blade_grid[k,i+1]).all(1))[0]       ## LL (2)
            p3_index = np.flatnonzero((blade == blade_grid[k+1,i+1]).all(1))[0]     ## LR (3)
            p4_index = np.flatnonzero((blade == blade_grid[k+1,i]).all(1))[0]       ## UR (4)

            ### ----------------------------    DUST    ---------------------------- #####
            ## If we close the TE in DUST, LAST POINT SHOULD NOT BE REPEATED !!!
            ## (For the last element in each span, connectivity will be linked to 1st element)

            # Connectivity will not change for the intermediate points
            if close_TE != 'yes' or (i_DUST+1)%n_chord_DUST!=0:

                p1_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k,i_DUST]).all(1))[0]         ## UL (1)   
                p2_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k,i_DUST+1]).all(1))[0]       ## LL (2)
                p3_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k+1,i_DUST+1]).all(1))[0]     ## LR (3)
                p4_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k+1,i_DUST]).all(1))[0]       ## UR (4)

            # Connectivity will be different for the last TE point !!!!!! 
            else:
                p1_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k,i_DUST]).all(1))[0]      ## UL (1)   
                p2_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k,0]).all(1))[0]    ## LL (2)
                p3_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k+1,0]).all(1))[0]         ## LR (3)
                p4_index_DUST = np.flatnonzero((blade_DUST == blade_DUST_grid[k+1,i_DUST]).all(1))[0]           ## UR (4)
        
            cell_connect = np.array([[p1_index, p2_index, p3_index, p4_index]])
            cell_connect_DUST = np.array([[p1_index_DUST, p2_index_DUST, p3_index_DUST, p4_index_DUST]])
        
            connectivity.append(cell_connect)
            connectivity_DUST.append(cell_connect_DUST)
            
    points = blade[:,1:]
    points_DUST = blade_DUST[:,1:]
    
    connectivity = np.concatenate([cell for cell in connectivity], axis=0)
    connectivity_DUST = np.concatenate([cell for cell in connectivity_DUST], axis=0)
    
    """ ##############         CREATE MESH       #############"""
    # Generate the DUST mesh #    
    cells = [('quad', connectivity_DUST)]
    mesh_DUST= meshio.Mesh(points_DUST, cells)
    
    # Generate the NVLM mesh for inspection #
    cells = [('quad', connectivity)]
    mesh= meshio.Mesh(points, cells)
    
    return mesh, mesh_DUST, connectivity_DUST, points_DUST
