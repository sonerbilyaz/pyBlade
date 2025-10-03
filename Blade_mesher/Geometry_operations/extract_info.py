import cadquery as cq

import numpy as np

def get_TE(edge_objects):
    """
        Gets the TE edge object from all of the edges
        
        Parameters
        ----------
        edge_objects    : List
            All edge objects of a section (includes TE)

        Returns
        -------
        TE              : Edge object
            TE edge object
        
        TE_upper_edge   : Edge object
            The edge at the upper surface adjacent to the TE 

        TE_lower_edge   : Edge object
            The edge at the lower surface adjacent to the TE 

    """

    # Get the midpoints of the edges
    edges_midpoints = np.concatenate([np.array(edge.positionAt(0.5).toTuple()).reshape((1,3)) for edge in edge_objects], axis=0)
                
    ## Select the TE whose midpoint has the highest x-coordinate ##
    TE_id = np.argmax(edges_midpoints[:,0])
    TE = edge_objects[TE_id]
    
    ## Get the upper and lower vertices of the TE
    TE_vertices = np.asarray([np.array(v.toTuple()) for v in TE.Vertices()])
    
    TE_v_up = TE_vertices[np.argmin(TE_vertices[:,1])]      ## Lowest y-coordinate will be upper TE vertex  
    TE_v_low = TE_vertices[np.argmax(TE_vertices[:,1])]     ## Highest y-coordinate will be upper TE vertex  

    ### Get the edges adjacent to the TE at the upper and lower surfaces
    for edge in edge_objects:
        v_array = np.asarray([np.array(v.toTuple()) for v in edge.Vertices()])

        # Get the upper edge
        if edge != TE and np.all(TE_v_up == v_array, axis = 1).any():
            TE_upper_edge = edge  
        # Get the lower edge
        elif edge != TE and np.all(TE_v_low == v_array, axis = 1).any():
            TE_lower_edge = edge  

    return TE, TE_upper_edge, TE_lower_edge

def extract_edges_open_TE(edges, TE_upper_edge, TE_lower_edge):
    """
        Extracts the upper and lower surfaces, seperately 
        
        Parameters
        ----------
        edges : List
            All edge objects of a section (WITHOUT TE!!)

        TE_upper_edge   : Edge object
            The edge at the upper surface adjacent to the TE 

        TE_lower_edge   : Edge object
            The edge at the lower surface adjacent to the TE 
            
        Returns
        -------
        lower_surface   : Wire object
            Lower surface as a single wire object

        upper_surface   : Wire object
            Upper surface as a single wire object
        
        LE              : array(3,)
            LE vertex coordinates
    """
    
    ### Find the LE vertex ###
    vertices = [np.array(v.toTuple()) for edge in edges for v in edge.vertices()]
    vertices = np.asarray(vertices)
    ## Lowest x-coordinate should be LE vertex
    LE = vertices[np.argmin(vertices[:,0])]
    
    for i, edge in enumerate(edges):
        if edge == TE_lower_edge:
            TE_low_id = i
     
    ######### IF THERE ARE ONLY UPPER AND LOWER SURFACE EDGES AFTER REMOVING THE TE #########
    if len(edges) == 2:
        upper_surface = TE_upper_edge
        lower_surface = TE_lower_edge

    ######### IF THERE ARE MORE EDGES AFTER REMOVING THE TE #########
    elif len(edges) > 2:       
        
        """
        The edge order is going CCW (seems so). Aim is to roll the elements such that
        the last element will be "TE_lower_edge", starting from "TE_upper_edge".
        """
        
        roll = len(edges) - TE_low_id % len(edges) -1
        ## Roll the edges list
        sorted_edges = edges[-roll:] + edges[:-roll]

        ## Make seperate wires for upper and lower surfaces from the sorted edges
        for i, edge in enumerate(sorted_edges):
            v_array = np.asarray([np.array(v.toTuple()) for v in edge.Vertices()])
            ## if LE vertex is reached, construct upper & lower surfaces
            if np.all(LE == v_array, axis = 1).any():  
                upper_surface = cq.Wire.assembleEdges(sorted_edges[:i+1])
                lower_surface = cq.Wire.assembleEdges(sorted_edges[i+1:])
                break
                        
    return upper_surface, lower_surface, LE
    
    