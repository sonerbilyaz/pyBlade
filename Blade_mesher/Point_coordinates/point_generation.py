import cadquery as cq # type: ignore
import numpy as np

from Distributions.distributions import airfoil_distribution

from Geometry_operations import modify, extract_info
from Geometry_operations.Rotation import Rotate

## Get coordinates ##
def get_coords(config, z_planes, collective_pitch, twist_local, chord_scale):
    
    ### Generate parametric points for airfoil section ###
    parametric_points, Node_ID_one_surf = airfoil_distribution(config["N_chord"], config["dist_section"])
    
    # Load the .stp file
    blade = cq.importers.importStep(f'{config["input_dir"]}/{config["file"]}')
    
    ALL_sections = []
    ALL_cross_sections = []
    ALL_sections_pitch = []
    ALL_sections_chord = []

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
        if config["CCW"] is False:
            edge_objects = [edge.mirror('ZY') for edge in edge_objects]
            edge_objects = edge_objects[::-1]
        
        ### Extract the Upper and Lower Surfaces Seperately
        upper_surface, lower_surface = extract_info.extract_surfaces(config, edge_objects)
        
        ## Assign the parametric points between 0 and 1 to the upper and lower surfaces
        parametric_points_up = parametric_points 
        Node_ID_up = Node_ID_one_surf
        
        parametric_points_low = parametric_points
        Node_ID_low = np.linspace(len(Node_ID_up)+1, 2*len(Node_ID_up)-1, len(Node_ID_up)-1)
        
        # Make sure it generates points starting from LE vertex  ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
        ## If the x-coord of the vertex is bigger than the other edge, reverse the order 
        reverse_upper_surf = upper_surface.positionAt(0).toTuple()[0] > upper_surface.positionAt(1).toTuple()[0]
        reverse_lower_surf = lower_surface.positionAt(0).toTuple()[0] > lower_surface.positionAt(1).toTuple()[0]

        if reverse_upper_surf and config["dist_section"] == 'cosine_LE':
            parametric_points_up = 1 - parametric_points_up
            
        if reverse_lower_surf and config["dist_section"] == 'cosine_LE':
            parametric_points_low = 1- parametric_points_low
            
        # Generate interpolated points along the upper and lower surfaces
        upper_points = np.asarray([np.array(upper_surface.positionAt(t).toTuple()) for t in parametric_points_up])
        lower_points = np.asarray([np.array(lower_surface.positionAt(t).toTuple()) for t in parametric_points_low])
        
        ## Increase pitch, if there are any ###
        upper_points = modify.pitch_increase(upper_points, config["initial_coll_pitch"] + collective_pitch)
        lower_points = modify.pitch_increase(lower_points, config["initial_coll_pitch"] + collective_pitch)

        # Insert Node ID  (From upper TE to lower TE)
        data_up  = np.insert(upper_points[:,:], 0, Node_ID_up, axis=1)    
        data_low = np.insert(lower_points[1:,:], 0, Node_ID_low,axis=1)
        
        data = np.append(data_up, data_low, axis=0)
        
        # Sort wrt ascending Node ID
        sorted_indices = np.argsort(data[:,0])
        data = data[sorted_indices]
        
        ## Make sectional changes, if there are any ###
        if np.any(twist_local!=0) or chord_scale!=1:
            data = modify.change_span(data, Node_ID_one_surf, twist, chord_scale)

        data_DUST = data.copy()

        # Check for closing the TE[:,1:]*1e-03
        if config["close_TE"] is True:
            data, pitch, chord = modify.close_TE_gap(data, Node_ID_one_surf)
            ## Coordinates will be different for DUST basic mesh !! 
            # Last repeated point at the TE should be removed 
            data_DUST = data[1:,:].copy()
            # After closing the TE, make sure that upper TE node and lower TE point is the same for the NVLM solver (1st and last point)
            data[-1,1:] = data[0,1:]
            
        
        # Convert from mm to meter
        data[:,1:] = data[:,1:]*1e-03
        data_DUST[:,1:] = data_DUST[:,1:]*1e-03
        
        ###     CHECK CCW !!!   ###
        # If CW ==> Reverse the sign of x-coordinates 
        # (mirror back to original)
        if config["CCW"] is False:
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
        ALL_sections_pitch.append(pitch)
        ALL_sections_chord.append(chord)

        data_DUST[:,0] = data_DUST[:,0][::-1]     ## Reverse node id
        ALL_sections_DUST.append(data_DUST[::-1]) ## Reverse points order
        
        # Cross sections append (Optional)
        ALL_cross_sections.append(cross_section.val())    
        
    ## Export all cross sections as a step file (Optional)
    all_sections_compound = cq.Compound.makeCompound(ALL_cross_sections)      
    
    return ALL_sections, ALL_sections_DUST, all_sections_compound, np.asarray(ALL_sections_pitch), np.asarray(ALL_sections_chord)