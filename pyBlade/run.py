import numpy as np

from .in_out.export_mesh import export

from .Distributions.distributions import spanwise_disribution as spanwise_planes
from .Point_coordinates import point_generation as points
from .Mesher.mesh import generate_mesh


def run(config, z_planes):
    
    #########   GENERATE POINTS  ##########
    ALL_sections, ALL_sections_DUST, all_sections_compound, pitch, chord = points.get_coords(config, z_planes, collective_pitch=0, twist_local=np.zeros_like(z_planes), chord_scale=1)
        
    #########   GENERATE MESH  #######
    mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, config["PANELS"]["close_TE"])

    #########   EXPORT   #######
    export(config, ALL_sections, z_planes, pitch, chord, mesh, mesh_DUST, connectivity_DUST, coordinates_DUST)

# def run_modify_section(config, z_planes):

    # for collective_pitch in config["coll_pitch_increment"]:
    #     for k in range(len(config["twist_root"])):
    #         for chord_scale in config["chord_scale"]:

    #             case_name = '{}_pitch_{:.2f}_twist_{:.2f}-{:.2f}_chord_{:.2f}'.format(case_prefix, float(collective_pitch), config["twist_root"][k] ,config["twist_tip"][k], chord_scale)
    #             dir_case = f'{config["output_dir"]}/{case_name}'

    #             ## If folder does not exist, create it
    #             if os.path.isdir(dir_case) is False:
    #                 os.mkdir(dir_case)


    #             twist_local = np.linspace(config["twist_root"][k], config["twist_tip"][k], len(z_planes))

    #             ###### GENERATE POINTS ######
    #             ALL_sections, ALL_sections_DUST, all_sections_compound, pitch, chord = points.get_coords(config, z_planes, collective_pitch, twist_local, chord_scale)
                    
    #             ###### GENERATE MESH AND EXPORT #######
    #             mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, n_blades, config["close_TE"])

    #             meshio.write(f'{dir_case}/{case_name}.vtk', mesh, file_format='vtk')

    #             #########    DUST EXPORT (for Basic Mesh)   #########
    #             DUST_dir = f'{dir_case}/DUST_output'

    #             if os.path.isdir(DUST_dir) is False:
    #                 os.mkdir(f'{DUST_dir}')

    #             meshio.write(f'{DUST_dir}/mesh_DUST.vtk', mesh_DUST, file_format='vtk')

    #             with open(f'{DUST_dir}/rr.dat', 'w') as file:
    #                 np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

    #             ## Export the DUST connectivity by switching from python index to dust index
    #             connectivity_DUST = connectivity_DUST + np.ones(connectivity_DUST.shape)

    #             ## !!!!!!! Connectivity will be reversed for the CW rotation !!!!!!!
    #             if config["CCW"] is False:
    #                 connectivity_DUST = connectivity_DUST[:,::-1]

    #             with open(f'{DUST_dir}/ee.dat', 'w') as file:
    #                 np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')

    #             ##################################   NVLM EXPORT     ##################################
    #             ALL_sections_array = np.concatenate(ALL_sections, axis=0)

    #             with open(f'{dir_case}/{case_name}.pts', 'w') as file:
    #                 file.write('######## Panel parameters ########\n')
    #                 file.write(f'type={surf_type}\n')
    #                 file.write(f'n_blades={n_blades}\n\n')
    #                 file.write(f'rotation_center={rotation_center}\n')
    #                 file.write(f'rotation_axis={rotation_axis}\n')
                    
    #                 file.write('n_span_all= {}\n'.format(len(ALL_sections)))
    #                 file.write('n_points={}\n'.format(len(ALL_sections[0][:,0])))
    #                 file.write('######## End of parameters ########\n')

    #                 ## Write the section coordinates ##                    
    #                 np.savetxt(file, ALL_sections_array, delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

###############################################################################################################################

