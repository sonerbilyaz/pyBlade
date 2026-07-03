import cadquery as cq # type: ignore
import numpy as np

from ..Distributions.distributions import airfoil_distribution

from ..Geometry_operations import extract, modify
from ..Geometry_operations.Rotation import Rotate

##  Create a function which creates a list of cq.Vectors from an array
def to_vec_list(coords_sec):
    """
        Converts a airfoil cross-sectional coordinates array into a list of cq.vectors
        to build a spline edge with "cq.Edge.makeSpline()"
        
        Parameters
        ----------
        coords_sec  : array (N,3) 
            Coordinates of the cross section 
            ( Goes from lower TE to upper TE, midpoint is the LE vertex !!! )

        Returns
        -------
        vectors     : list (N)
            List of cq.vectors

    """
    return [cq.Vector(float(x), float(y), float(z)) for x, y, z in coords_sec]


def create_wire(coords_sec):
    """
        Creates a wire object composed of upper surface and lower surface edges

        Parameters
        ----------
        coords_sec  : array (N,3) 
            Coordinates of the cross section 
            ( Goes from lower TE to upper TE, midpoint is the LE vertex !!! )

        Returns
        -------
        wire        : cq.Wire object 
            Wire object consists of lower and upper surface edges
    """

    ## Seperate upper and lower surfaces
    LE_ID = coords_sec.shape[0] // 2

    # Create a list of vectors to make a spline edge for upper and lower surf
    coords_sec_up  = to_vec_list(coords_sec[LE_ID:])    # cq.Vector
    coords_sec_low = to_vec_list(coords_sec[:LE_ID+1])  # cq.Vector
    
    # Convert them into edges with spline operation
    upper_edge = cq.Edge.makeSpline(coords_sec_up)
    lower_edge = cq.Edge.makeSpline(coords_sec_low)

    ##  Create a single wire section ## 
    wire = cq.Wire.assembleEdges([lower_edge, upper_edge])

    return wire


def create_Blade(config, sections, collective_pitch, twist_local, chord_scale):
    """
        Creates the 3D blade geometry from the blade elements 
        with Loft operation

        Parameters
        ----------
        config              : dict
            Configuration parameters

        sections            : array (N,)
            Radial locations of the sections

        collective_pitch    : list 
            Collective pitch increment (deg)

        twist_local         : array (N,)
            Individual twist of the sections (deg). Generally zeros for uniform change along span

        chord_scale         : list
            Cross section scale

        Returns
        -------
        blade               : cq.occ_impl.shapes.Solid object
            Solid lofted shape of the blade
    """
    
    ### Generate parametric points for airfoil section ###
    parametric_points, Node_ID_one_surf = airfoil_distribution(int(config["PANEL"]["N_chord"]), config["PANEL"]["dist_chord"])
    
    # Load the .stp file
    blade = cq.importers.importStep(f'{config["FILE"]["stp_file"]}')

    ALL_sections_wire = []
    ALL_sections_pitch = []
    ALL_sections_chord = []

    for sec, twist in zip(sections, twist_local):
        # Create a section of the blade using the defined z-coord
        cross_section=blade.section(height=sec)
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
        upper_surface, lower_surface = extract.get_surfaces(config, edge_objects)
        
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
        upper_points = modify.pitch_increase(upper_points, collective_pitch)    # (N,3)
        lower_points = modify.pitch_increase(lower_points, collective_pitch)    # (M,3)

        # Insert Node ID  (From upper TE to lower TE)
        coords_up  = np.insert(upper_points[:,:], 0, Node_ID_up, axis=1)        # (N,4)
        coords_low = np.insert(lower_points[1:,:], 0, Node_ID_low,axis=1)       # (M,4)

        coords_sec = np.append(coords_up, coords_low, axis=0)                   # (S=M+N,4)
        
        # Sort wrt ascending Node ID
        sorted_indices = np.argsort(coords_sec[:,0])                            # (S,4)
        coords_sec = coords_sec[sorted_indices]                                 # (S,4)
        
        ## Make sectional changes, if there are any ###
        if np.any(twist_local!=0) or np.any(chord_scale!=1):
            coords_sec = modify.change_span(coords_sec, Node_ID_one_surf, twist, chord_scale)

        #####    Check for closing the TE   #####
        if config["SURFACE"]["close_TE"] in ['yes', True, 'Yes']:
            coords_sec, pitch, chord = modify.close_TE_gap(coords_sec, Node_ID_one_surf)
            # After closing the TE, make sure that upper TE node and lower TE point is the same for the NVLM solver (1st and last point)
            coords_sec[-1,1:] = coords_sec[0,1:]
                
        ##########     CHECK CCW !!!   ##########
        # If CW ==> Reverse the sign of x-coordinates (mirror back to original)
        if config["SURFACE"]["rotation"] == "CW":
            coords_sec[:,1] = -coords_sec[:,1]

        ##############      Wire Generation        ##############
        ##  First, remove the Node_ID    (Remove 1st column) ##
        coords_sec = np.delete(coords_sec, np.s_[0], 1)
        ## Create the wire ##
        wire_sec = create_wire(coords_sec)

        ## Append the data ##
        ALL_sections_wire.append(wire_sec)
        ALL_sections_pitch.append(pitch)
        ALL_sections_chord.append(chord)

    ## Create the blade (Loft) ##
    blade = cq.Solid.makeLoft(ALL_sections_wire, ruled = True)

    return blade