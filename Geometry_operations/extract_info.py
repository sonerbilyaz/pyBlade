import cadquery as cq
import numpy as np

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
    
    return upper_surface, lower_surface, LE_coordinates