import numpy as np, ast, cadquery as cq         # type: ignore

from .in_out import export

from .Distributions.distributions import spanwise_disribution as spanwise_planes
from .Point_coordinates import point_generation as points
from .Mesher.mesh import generate_mesh
from .Blade_Generate.blade import create_Blade
from .Blade_Generate import mesh_blade

def run_mesh(config):
    
    # Load the .stp file
    blade = cq.importers.importStep(f'{config["FILE"]["stp_file"]}')
    
    ## Create spanwise planes ###
    z_planes = spanwise_planes(config["PANEL"]["z_min"], config["PANEL"]["z_max"], config["PANEL"]["N_span"], config["PANEL"]["dist_span"], config["PANEL"]["r_R"])

    collect_pitch = np.array(ast.literal_eval(config["PANEL"]["collect_pitch"])).astype(float)
    scale = np.array(ast.literal_eval(config["PANEL"]["scale"])).astype(float)
    sweep_angle = np.array(ast.literal_eval(config["PANEL"]["sweep"])).astype(float)
    dihedral_angle = np.array(ast.literal_eval(config["PANEL"]["dihedral"])).astype(float)
    twist_tip = np.array(ast.literal_eval(config["PANEL"]["twist_tip"])).astype(float)
    twist_root = np.array(ast.literal_eval(config["PANEL"]["twist_root"])).astype(float)

    twist_span = np.concatenate([[twist_tip], [twist_root]], axis=0)

    for coll_pitch in collect_pitch:
        for chord_scale in scale:
            for sweep in sweep_angle:
                for dihedral in dihedral_angle:
                    for i in range(twist_span.shape[1]):
                        twist_tip, twist_root = twist_span[0,i], twist_span[1,i]
                        
                        no_change_twist = twist_root == 0 and twist_tip == 0
                        no_change_uniform = config["PANEL"]["collect_pitch"] == '[0]' and config["PANEL"]["scale"] == '[1]' and config["PANEL"]["sweep"] == '[0]' and config["PANEL"]["dihedral"] == '[0]'
                        
                        if no_change_twist and no_change_uniform:
                            case_name = config["FILE"]["name_tag"]                
                        else:
                            case_name = '{}_{:.1f}deg_coll_pitch_{:.2f}_scale_{:.1f}deg_sweep_{:.1f}deg_dihedr_{}-{}deg_twist'.format(config["FILE"]["name_tag"], coll_pitch, chord_scale, sweep, dihedral, twist_tip, twist_root)

                        ######  GENERATE POINTS  ######
                        ALL_sections, ALL_sections_DUST, pitch, chord = points.get_coords(blade, config, z_planes, coll_pitch, chord_scale, sweep, dihedral, twist_tip, twist_root)
                            
                        ######  GENERATE MESH   #######
                        mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, config["SURFACE"]["close_TE"])

                        #########   EXPORT   #######
                        export.export_mesh(config, case_name, ALL_sections, z_planes, pitch, chord, mesh, mesh_DUST, connectivity_DUST, coordinates_DUST)

def run_blade(config):

    ## First, create sections ##
    sections = spanwise_planes(config["GENERATE_SURFACE"]["z_min_sec"], config["GENERATE_SURFACE"]["z_max_sec"], config["GENERATE_SURFACE"]["n_sec"], config["GENERATE_SURFACE"]["dist_sec"], config["GENERATE_SURFACE"]["r_R_sec"])
    
    collect_pitch_sec = np.array(ast.literal_eval(config["GENERATE_SURFACE"]["collect_pitch_sec"])).astype(float)
    scale_sec = np.array(ast.literal_eval(config["GENERATE_SURFACE"]["scale_sec"])).astype(float)

    ##  Create spanwise planes (FOR MESHING)  ###
    z_planes = spanwise_planes(config["PANEL"]["z_min"], config["PANEL"]["z_max"], config["PANEL"]["N_span"], config["PANEL"]["dist_span"], config["PANEL"]["r_R"])

    for coll_pitch_sec in collect_pitch_sec:
        for chord_scale_sec in scale_sec:

            if config["GENERATE_SURFACE"]["collect_pitch_sec"] == '[0]' and config["GENERATE_SURFACE"]["scale_sec"] == '[1]':
                case_name_sec = '{}_{}n_sec'.format(config["FILE"]["name_tag"], config["GENERATE_SURFACE"]["n_sec"])
            else:
                case_name_sec = '{}_{}n_sec_{:.2f}_delta_coll_pitch_{:.2f}_scale'.format(config["FILE"]["name_tag"], config["GENERATE_SURFACE"]["n_sec"], coll_pitch_sec, chord_scale_sec)

            ######  GENERATE BLADE  ######
            blade = create_Blade(config, sections, coll_pitch_sec, np.zeros_like(sections), scale_sec)

            #####     If mesh will be generated, DO THAT  #####
            if config["GENERATE_SURFACE"]["generate_mesh"] in ['yes', True, 'Yes']:
                
                ### Convert the blade into a Workplane object to mesh it
                blade_WP = cq.Workplane(blade)      # cq.Workplane

                ######  GENERATE POINTS  (LE is already defined, and TE is already closed !!)######
                ALL_sections, ALL_sections_DUST = mesh_blade.mesh(blade_WP, config, z_planes, coll_pitch_sec, np.zeros_like(z_planes), scale_sec)
                
                ######  GENERATE MESH   #######
                mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, close_TE='yes')

                #########   EXPORT MESH   #######
                export.export_blade_mesh(config, case_name_sec, ALL_sections, mesh, mesh_DUST, connectivity_DUST, coordinates_DUST)

            #########   EXPORT BLADE   #######
            export.export_blade(config, case_name_sec, blade)

