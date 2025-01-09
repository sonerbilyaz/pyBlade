import cadquery as cq
import numpy as np

from mesh_from_CAD.Distributions.distributions import airfoil_distribution

from mesh_from_CAD.Geometry_operations import modify as modify
from mesh_from_CAD.Geometry_operations import extract_info as identify

## Get coordinates ##
def get_coords(stp_file, num_points, dist_airfoil, z_planes, remove_TE, close_TE):    
    
    ### Generate parametric points for airfoil section ###
    parametric_points, Node_ID_one_surf = airfoil_distribution(num_points, dist_airfoil)
    
    # Load the .stp file
    blade = cq.importers.importStep(stp_file)
    
    ALL_sections = []
    ALL_cross_sections = []
    
    for z in z_planes:
        # Create a section of the blade using the defined z-coord
        cross_section=blade.section(height=z)
        # Get all edges from the cross section
        edges = cross_section.edges()
        # Access the edge objects from the Workplane
        edge_objects = edges.objects
        
        # Should we remove TE ??
        if remove_TE is True:
            filtered_edges, TE_upper, TE_lower = modify.remove_TE(edge_objects)
        
        if remove_TE is False:
            _, TE_upper, TE_lower  = modify.remove_TE(edge_objects)
            filtered_edges = edge_objects
        
        # Identify the edges
        upper_surface, lower_surface, LE_coords = identify.extract_edges(filtered_edges, TE_upper, TE_lower)
        
        ## Assign the parametric points between 0 and 1 to the upper and lower surfaces
        parametric_points_up = parametric_points 
        Node_ID_up = Node_ID_one_surf
        
        parametric_points_low = parametric_points
        Node_ID_low = np.linspace(len(Node_ID_up)+1, 2*len(Node_ID_up)-1, len(Node_ID_up)-1)
        
        # Check it starts generating points from the starting vertex (LE)
        if np.array_equal(np.round(upper_surface.positionAt(0).toTuple(),decimals=2).reshape(1,3), np.round(LE_coords,decimals=2)) is False and dist_airfoil == 'cosine_LE':
            # If the orientation is wrong, reverse it to correct ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
            parametric_points_up = 1 - parametric_points_up
            # 
            # print(f'upper surface parametric is reversed at z={z}mm. LE = {upper_surface.positionAt(parametric_points_up[0]).toTuple()}')
            
        if np.array_equal(np.round(lower_surface.positionAt(0).toTuple(),decimals=2).reshape(1,3), np.round(LE_coords,decimals=2)) is False and dist_airfoil == 'cosine_LE':
            # If the orientation is wrong, reverse it to correct ==> !! ONLY REQUIRED FOR COSINE_LE DISTRIBUTION !!
            parametric_points_low = 1- parametric_points_low
            # 
            # print(f'lower surface parametric is reversed at z={z}mm LE = {lower_surface.positionAt(parametric_points_low[0]).toTuple()}')
            
        # Generate interpolated points along the upper and lower surfaces
        upper_points, lower_points = [], []
        for t in parametric_points_up:
            # Place points on each edge using normalized parameter t
            upper_point  = upper_surface.positionAt(t).toTuple()
            upper_points.append(upper_point)
                
        for k in parametric_points_low:
            lower_point  = lower_surface.positionAt(k).toTuple()
            lower_points.append(lower_point)
            
        # Create coords
        upper_points = np.array(upper_points)
        lower_points = np.array(lower_points)
        
        # Insert Node ID  (From upper TE to lower TE)
        data_up  = np.insert(upper_points[:,:], 0, Node_ID_up, axis=1)    
        data_low = np.insert(lower_points[1:,:], 0, Node_ID_low,axis=1)
        
        data = np.append(data_up, data_low, axis=0)
        
        # Sort wrt ascending Node ID
        sorted_indices = np.argsort(data[:,0])
        data = data[sorted_indices]
        
        # Check for closing the TE[:,1:]*1e-03
        if close_TE is True:
            data = modify.close_TE_gap(data, Node_ID_one_surf, n=10)
            # After closing the TE, remove extra TE point, which is the last element
            data = np.delete(data, -1, axis=0)
        
        # Convert from mm to meter
        data[:,1:] = data[:,1:]*1e-03
        
        ### Append ALL data ###
        ALL_sections.append(data)
        # Cross sections append (Optional)
        ALL_cross_sections.append(cross_section.val())    
        
    ## Export all cross sections as a step file (Optional)
    all_sections_compound = cq.Compound.makeCompound(ALL_cross_sections)      
    
    return ALL_sections, all_sections_compound