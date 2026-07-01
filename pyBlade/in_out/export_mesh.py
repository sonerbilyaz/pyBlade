import meshio
import os, numpy as np

def export(config, ALL_sections: list, z_planes, pitch, chord, mesh, mesh_DUST, connectivity_DUST, coordinates_DUST, DUST_export = True):

    # Create output directory if it is absent #
    if os.path.isdir(config["FILE"]["output_dir"]) is False:
        os.mkdir(config["FILE"]["output_dir"])

    ## If case folder does not exist, create it
    case_dir = f'{config["FILE"]["output_dir"]}/{config["FILE"]["name_tag"]}'
    if os.path.isdir(case_dir) is False:
        os.mkdir(case_dir)

    ###################     Normal mesh export     ###################
    meshio.write(f'{case_dir}/{config["FILE"]["name_tag"]}.vtk', mesh, file_format='vtk')

    ## Generate array for ALL_sections (NODE_ID + coordinates)
    ALL_sections_array = np.concatenate(ALL_sections, axis=0)

    with open(f'{case_dir}/{config["FILE"]["name_tag"]}.pts', 'w') as file:
        file.write('######## Panel parameters ########\n')
        # file.write(f'type={surf_type}\n')
        file.write(f'n_blades={int(config["SURFACE"]["n_blades"])}\n\n')
        # file.write(f'rotation_center={rotation_center}\n')
        # file.write(f'rotation_axis={rotation_axis}\n')
        
        file.write('n_span_all= {}\n'.format(len(ALL_sections)))
        file.write('n_points={}\n'.format(len(ALL_sections[0][:,0])))
        file.write('######## End of parameters ########\n')
        
        ## Write the section coordinates ##                    
        np.savetxt(file, ALL_sections_array, delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

    ###    Export pitch and chord   ###
    np.savetxt(f'{case_dir}/{config["FILE"]["name_tag"]}_pitch_chord.txt', 
               np.concatenate([z_planes[:, None]/float(config["PANEL"]["z_max"]), chord[:, None], pitch[:,None]], axis=1), header='r_R' + '\t' + 'chord (mm)' + '\t' +'pitch(deg)', delimiter='\t')
    
    if DUST_export:

        ###################       DUST export     ###################
        DUST_dir = f'{case_dir}/DUST_output'
        if os.path.isdir(DUST_dir) is False:
            os.mkdir(f'{DUST_dir}')

        ## vtk export for visualisation
        meshio.write(f'{DUST_dir}/{config["FILE"]["name_tag"]}_DUST.vtk', mesh_DUST, file_format='vtk')     

        ##  DUST coordinates  ##
        with open(f'{DUST_dir}/rr.dat', 'w') as file:
            np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

        ## Export the DUST connectivity by switching from python index to dust index
        connectivity_DUST = connectivity_DUST + np.ones(connectivity_DUST.shape)
        ## !!! Connectivity will be reversed for the CW rotation !!!
        if config["IDENTIFY"]["rotation"] == "CW":
            connectivity_DUST = connectivity_DUST[:,::-1]
        
        ##  DUST connectivity  ##
        with open(f'{DUST_dir}/ee.dat', 'w') as file:
            np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')

    # ## Export all cross sections (Optional) ##
    # cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")


def export_modified(config, case_name, ALL_sections: list, z_planes, pitch, chord, mesh, mesh_DUST, connectivity_DUST, coordinates_DUST, DUST_export = True):

    # Create output directory if it is absent #
    if os.path.isdir(config["FILE"]["output_dir"]) is False:
        os.mkdir(config["FILE"]["output_dir"])

    ## If case folder does not exist, create it
    case_dir = f'{config["FILE"]["output_dir"]}/{case_name}'
    if os.path.isdir(case_dir) is False:
        os.mkdir(case_dir)

    ###################     Normal mesh export     ###################
    meshio.write(f'{case_dir}/{case_name}.vtk', mesh, file_format='vtk')

    ## Generate array for ALL_sections (NODE_ID + coordinates)
    ALL_sections_array = np.concatenate(ALL_sections, axis=0)

    with open(f'{case_dir}/{config["FILE"]["name_tag"]}.pts', 'w') as file:
        file.write('######## Panel parameters ########\n')
        # file.write(f'type={surf_type}\n')
        file.write(f'n_blades={int(config["SURFACE"]["n_blades"])}\n\n')
        # file.write(f'rotation_center={rotation_center}\n')
        # file.write(f'rotation_axis={rotation_axis}\n')
        
        file.write('n_span_all= {}\n'.format(len(ALL_sections)))
        file.write('n_points={}\n'.format(len(ALL_sections[0][:,0])))
        file.write('######## End of parameters ########\n')
        
        ## Write the section coordinates ##                    
        np.savetxt(file, ALL_sections_array, delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

    ###    Export pitch and chord   ###
    np.savetxt(f'{case_dir}/{config["FILE"]["name_tag"]}_pitch_chord.txt', 
               np.concatenate([z_planes[:, None]/float(config["PANEL"]["z_max"]), chord[:, None], pitch[:,None]], axis=1), header='r_R' + '\t' + 'chord (mm)' + '\t' +'pitch(deg)', delimiter='\t')

    if DUST_export:

        ###################       DUST export     ###################
        DUST_dir = f'{case_dir}/DUST_output'
        if os.path.isdir(DUST_dir) is False:
            os.mkdir(f'{DUST_dir}')

        ## vtk export for visualisation
        meshio.write(f'{DUST_dir}/{config["FILE"]["name_tag"]}_DUST.vtk', mesh_DUST, file_format='vtk')

        ##  DUST coordinates  ##
        with open(f'{DUST_dir}/rr.dat', 'w') as file:
            np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

        ## Export the DUST connectivity by switching from python index to dust index
        connectivity_DUST = connectivity_DUST + np.ones(connectivity_DUST.shape)
        ## !!! Connectivity will be reversed for the CW rotation !!!
        if config["IDENTIFY"]["rotation"] == "CW":
            connectivity_DUST = connectivity_DUST[:,::-1]
        
        ##  DUST connectivity  ##
        with open(f'{DUST_dir}/ee.dat', 'w') as file:
            np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')

    # ## Export all cross sections (Optional) ##
    # cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")
