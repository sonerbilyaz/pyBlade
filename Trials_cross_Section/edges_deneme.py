from IPython import get_ipython
get_ipython().run_line_magic('clear', '')

import numpy as np, os 

curr_dir = os.getcwd()

import cadquery as cq

os.chdir('../')
from Geometry_operations.modify import remove_TE
os.chdir(curr_dir)


stp_file = '../Test_Cases/12x6 ClarkY/12x6_ClarkY-1_Blade_LE_at_the_mid_face_little_modification.stp'
# import blade
blade = cq.importers.importStep(stp_file)

# Get section
z= 111           # sec location
z_max = 152     # max span location
section = blade.section(height=z)
cq.exporters.export(section, 'working section.step')
# Get edge objects
edges = section.edges().objects
# remove TE
filtered_edges, TE_upper, TE_lower = remove_TE(edges)

####### Identify LE #######
### Compare the midpoint coordinates of the each edge. Lowest x-coord midpoint should adjacent to the LE vertex ###
midpoint_edges = np.zeros((len(filtered_edges),3))
for i,edge in enumerate(filtered_edges):
    midpoint_edges[i] = np.array(edge.positionAt(0.5).toTuple()).reshape(1,3)
    
# Find the ID and edge which has the lowest x-coordinate among the midpoints. This ID is the edge id.
index = np.argmin(midpoint_edges[:,0])
LE_edge = filtered_edges[index]

""" IF LE IS NOT SPECIFIED AND IT IS SOMEWHERE AT THE MIDPOINT OF THE LEFTMOST EDGE """
## Get the LE coordinates which are at the midpoint of the leftmost edge
LE_coordinates = np.array(LE_edge.positionAt(0.5).toTuple()).reshape(1,3)
## Create the LE vertex for splitting
LE_vertex = cq.Vertex.makeVertex(LE_coordinates[0,0], LE_coordinates[0,1], LE_coordinates[0,2])

""" 
We are going to create a single combined wire which contains the ALL edges 
(upper & lower surface) except TE. To do this, we need to convert each edge into a seperate wire.
After doing that, we will split the wire at the LE vertex loc to create LE. By doing this,
we are creating the LE and sorting the new split edges so that endpoints will be the end edges

NOTE: If possible, stitching each adjacent wire will be perfect but couldnt do that. 
Command => cq.Wire.stitch(wire1, wire2) 
"""
# Convert each edge into a wire
filtered_wires = [cq.Wire.assembleEdges([edge]) for edge in filtered_edges]
# Merge the wires to get the overall cross section as a wire
combined_wire = cq.Wire.combine(filtered_wires, tol=1e-01)

# Split the whole wire at the LE vertex location
split = combined_wire[0].split(LE_vertex)
# Get the edges as a list. They are ordered now
split_edges = [edge for edge in split]

######## Create combined wires for UPPER SURFACE and LOWER SURFACE ###########
combined_wire_ALL = cq.Wire.combine(split_edges)[0]

## Find the LE vertex and its index on the new sorted combined wire
for i,vertex in enumerate(combined_wire_ALL.Vertices()):
    if np.array_equal(np.array(vertex.toTuple()).reshape(1,3), LE_coordinates):
        LE_vertex_index = i

## Combine the edges as a wire before and after the LE vertex ##
for i in range(len(combined_wire_ALL.Edges())):
    if i == LE_vertex_index-1:
        split_1 = cq.Wire.combine(combined_wire_ALL.Edges()[0:i+1])[0]
        split_2= cq.Wire.combine(combined_wire_ALL.Edges()[i+1:])[0]
        break

### IDENTIFY the UPPER SURFACE and LOWER SURFACE seperately ###
wires = [split_1, split_2]
for wire in wires:
    if np.array_equal(np.round(wire.positionAt(0).toTuple(), decimals=1), np.round(TE_lower, decimals=1)) or np.array_equal(np.round(wire.positionAt(1).toTuple(), decimals=1), np.round(TE_lower, decimals=1)):
        lower_surface= wire
    
    if np.array_equal(np.round(wire.positionAt(0).toTuple(),decimals=1), np.round(TE_upper, decimals=1)) or np.array_equal(np.round(wire.positionAt(1).toTuple(), decimals=1), np.round(TE_upper, decimals=1)):
        upper_surface= wire

print(f'LE {LE_coordinates}')
print(f'upper {upper_surface.positionAt(0).toTuple()}')
print(f'upper {upper_surface.positionAt(1).toTuple()}')
print(f'lower {lower_surface.positionAt(0).toTuple()}')
print(f'lower {lower_surface.positionAt(1).toTuple()}')
