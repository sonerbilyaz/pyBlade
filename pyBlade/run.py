import numpy as np, ast

from .in_out import export_mesh as export

from .Distributions.distributions import spanwise_disribution as spanwise_planes
from .Point_coordinates import point_generation as points
from .Mesher.mesh import generate_mesh


def run(config, z_planes):
    
    #########   GENERATE POINTS  ##########
    ALL_sections, ALL_sections_DUST, all_sections_compound, pitch, chord = points.get_coords(config, z_planes, collective_pitch=0, twist_local=np.zeros_like(z_planes), chord_scale=1)
        
    #########   GENERATE MESH  #######
    mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, config["PANEL"]["close_TE"])

    #########   EXPORT   #######
    export.export(config, ALL_sections, z_planes, pitch, chord, mesh, mesh_DUST, connectivity_DUST, coordinates_DUST)

def run_modify_planform(config, z_planes):

    collect_pitch = np.array(ast.literal_eval(config["MODIFY"]["collect_pitch"])).astype(float)
    scale = np.array(ast.literal_eval(config["MODIFY"]["scale"])).astype(float)

    for coll_pitch in collect_pitch:
        for chord_scale in scale:

            case_name = '{}_{:.2f}_delta_coll_pitch_{:.2f}_scale'.format(config["FILE"]["name_tag"], coll_pitch, chord_scale)

            ######  GENERATE POINTS  ######
            ALL_sections, ALL_sections_DUST, all_sections_compound, pitch, chord = points.get_coords(config, z_planes, coll_pitch, np.zeros_like(z_planes), chord_scale)
                
            ######  GENERATE MESH   #######
            mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, config["PANEL"]["close_TE"])

            #########   EXPORT   #######
            export.export_modified(config, case_name, ALL_sections, z_planes, pitch, chord, mesh, mesh_DUST, connectivity_DUST, coordinates_DUST)