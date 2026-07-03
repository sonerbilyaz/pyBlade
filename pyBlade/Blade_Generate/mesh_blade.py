import cadquery as cq # type: ignore
import numpy as np

from ..Distributions.distributions import airfoil_distribution

from ..Geometry_operations import extract, modify
from ..Geometry_operations.Rotation import Rotate

## Get coordinates ##
def mesh(blade, config, z_planes, collective_pitch, twist_local, chord_scale):
    
    ### Generate parametric points for airfoil section ###
    parametric_points, Node_ID_one_surf = airfoil_distribution(int(config["PANEL"]["N_chord"]), config["PANEL"]["dist_chord"])
    
    ALL_sections = []
    ALL_sections_DUST = []

    for z, twist in zip(z_planes, twist_local):
        # Create a section of the blade using the defined z-coord
        cross_section=blade.section(height=z)
        # Get all edges from the cross section
        edges = cross_section.edges()
        # Access the edge objects from the Workplane
        edge_objects = edges.objects

        ## First, check the rotation direction (CCW / CW ??)
        # CW ==> Reverse the sign of x-coord (mirror)
        if config["SURFACE"]["rotation"] == "CW":
            edge_objects = [edge.mirror('ZY') for edge in edge_objects]
            edge_objects = edge_objects[::-1]
        
        ### Extract the Upper and Lower Surfaces Seperately

        # Get the midpoints of the edges
        edges_midpoints = np.concatenate([np.array(edge.positionAt(0.5).toTuple()).reshape((1,3)) for edge in edge_objects], axis=0)

        ## Lower midpoint y-coordinate will be upper surface 
        if edges_midpoints[0,1] < edges_midpoints[1,1]:
            upper_surface = edge_objects[0]     ## Edge object
            lower_surface = edge_objects[1]     ## Edge object

        else:
            upper_surface = edge_objects[1]     ## Edge object
            lower_surface = edge_objects[0]     ## Edge object
        
        ## Assign the parametric points between 0 and 1 to the upper and lower surfaces
        parametric_points_up = parametric_points 
        Node_ID_up = Node_ID_one_surf
        
        parametric_points_low = parametric_points
        Node_ID_low = np.linspace(len(Node_ID_up)+1, 2*len(Node_ID_up)-1, len(Node_ID_up)-1)
        
        # Make sure it generates points starting from LE vertex  ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
        ## If the x-coord of the vertex is bigger than the other edge, reverse the order 
        reverse_upper_surf = upper_surface.positionAt(0).toTuple()[0] > upper_surface.positionAt(1).toTuple()[0]
        reverse_lower_surf = lower_surface.positionAt(0).toTuple()[0] > lower_surface.positionAt(1).toTuple()[0]

        if reverse_upper_surf and config["PANEL"]["dist_chord"] == 'cosine_LE':
            parametric_points_up = 1 - parametric_points_up
            
        if reverse_lower_surf and config["PANEL"]["dist_chord"] == 'cosine_LE':
            parametric_points_low = 1- parametric_points_low
            
        # Generate interpolated points along the upper and lower surfaces
        upper_points = np.asarray([np.array(upper_surface.positionAt(t).toTuple()) for t in parametric_points_up])
        lower_points = np.asarray([np.array(lower_surface.positionAt(t).toTuple()) for t in parametric_points_low])
        
        ## Increase pitch, if there are any ###
        upper_points = modify.pitch_increase(upper_points, collective_pitch)
        lower_points = modify.pitch_increase(lower_points, collective_pitch)

        # Insert Node ID  (From upper TE to lower TE)
        data_up  = np.insert(upper_points[:,:], 0, Node_ID_up, axis=1)    
        data_low = np.insert(lower_points[1:,:], 0, Node_ID_low,axis=1)
        
        data = np.append(data_up, data_low, axis=0)
        
        # Sort wrt ascending Node ID
        sorted_indices = np.argsort(data[:,0])
        data = data[sorted_indices]
        
        ## Make sectional changes, if there are any ###
        if np.any(twist_local!=0) or np.any(chord_scale!=1):
            data = modify.change_span(data, Node_ID_one_surf, twist, chord_scale)

        # For DUST output, th elast repeated point at the TE should be removed 
        data_DUST = data[1:,:].copy()

        # Convert from mm to meter
        data[:,1:] = data[:,1:]*1e-03
        data_DUST[:,1:] = data_DUST[:,1:]*1e-03
        
        ###     CHECK CCW !!!   ###
        # If CW ==> Reverse the sign of x-coordinates 
        # (mirror back to original)
        if config["SURFACE"]["rotation"] == "CW":
            data[:,1] = -data[:,1]
            data_DUST[:,1] = -data_DUST[:,1]

        ### ADJUST THE ORIENTATION !! (z-axis rotation same as DUST!!)
        data[:,1:] = Rotate(data[:,1:], -90, axis='x')
        data[:,1:] = Rotate(data[:,1:], -90, axis='z')
        data_DUST[:,1:] = Rotate(data_DUST[:,1:], -90, axis='x')
        data_DUST[:,1:] = Rotate(data_DUST[:,1:], -90, axis='z')
        
        ### Append ALL data ###
        ## Sectional points should go from lower to upper TE !!!
        data[:,0] = data[:,0][::-1]     ## Reverse node id
        ALL_sections.append(data[::-1]) ## Reverse points order

        data_DUST[:,0] = data_DUST[:,0][::-1]     ## Reverse node id
        ALL_sections_DUST.append(data_DUST[::-1]) ## Reverse points order
                    
    return ALL_sections, ALL_sections_DUST