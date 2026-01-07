from timeit import default_timer as timer
start_time = timer()  # Start the timer

import numpy as np, os, sys 

# Get the absolute path to the directory and add it to the sys path for relative imports
package_dir = os.path.abspath(os.path.dirname(__file__))
if package_dir not in sys.path:
    sys.path.append(package_dir)

from Distributions.distributions import spanwise_disribution as spanwise_planes
from Point_coordinates import point_generation as points
from Mesher.mesh import generate_mesh

"""################################ INPUTS #######################################"""
# File paths and the STEP file
working_dir = '../Runs/12x6_ClarkY'
stp_file = f'{working_dir}/12x6_ClarkY-1_Blade.stp'

output_dir = f'{working_dir}/output'
pts_filename = '12x6_ClarkY'

find_LE = True      ## Should we identify the LE ?? If so ==> True
###############     Panel Parameters     ###############
### Chordwise Distribution ###
N_chord = 23                    ## Upper and Lower surf separately!!
dist_airfoil = 'cosine_LE'

### Spanwise Distribution ###
N_span= 55
z_min, z_max = 23.5, 151.5          ## in mm

dist_spanwise = 'cosine_TIP'
r_R = 0.82                      ## Span location to start the cosine_TIP 

### TE Modification ###
remove_TE = True            # Should we remove TE ??
close_TE = True             # Should we close the TE gap ??

### Pitch Modification (deg) ###
collective_pitch_increment = [0,10,20]

### Spanwise local Changes ###
taper_local_root = [1,1.5]
taper_local_tip = [1,1.5]

twist_local_root = [0,10]
twist_local_tip = [0,0]

### .pts Inputs ###
n_blades = 1        # Number of blades
surf_type = 'propeller'     # Surface type (wing/propeller)
rotation_center = [0, 0, 0]
rotation_axis = [0, -1, 0]

"""###############################################################################"""
# Create output directory if it is absent #
if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

### Create flag for the output file name ###
cond_no_change = np.all(np.array(collective_pitch_increment) == 0) and np.all(np.array(taper_local_tip)==1) and np.all(np.array(taper_local_root)==1) and np.all(np.array(twist_local_root)==0) and np.all(np.array(twist_local_tip)==0)

### Create spanwise planes ###
z_planes = spanwise_planes(z_min, z_max, N_span, dist_spanwise, r_R)

if cond_no_change:
    ###### GENERATE POINTS ######
    ALL_sections, ALL_sections_DUST, all_sections_compound = points.get_coords(stp_file, N_chord, dist_airfoil, z_planes, close_TE, 0, twist_local=np.zeros_like(z_planes), taper_local=np.ones_like(z_planes), find_LE=find_LE)
        
    ###### GENERATE MESH AND EXPORT #######
    mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, n_blades, close_TE)

    import meshio
    meshio.write(f'{output_dir}/{pts_filename}.vtk', mesh, file_format='vtk')

    """####################  DUST EXPORT (for Basic Mesh)   ###########################"""
    DUST_dir = f'{working_dir}/DUST_output'
    if os.path.isdir(DUST_dir) is False:
        os.mkdir(f'{DUST_dir}')

    meshio.write(f'{DUST_dir}/mesh_DUST.vtk', mesh_DUST, file_format='vtk')

    with open(f'{DUST_dir}/rr.dat', 'w') as file:
        np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

    ## Export the DUST connectivity by switching from python index to dust index
    connectivity_DUST = connectivity_DUST + np.ones(connectivity_DUST.shape)
    with open(f'{DUST_dir}/ee.dat', 'w') as file:
        np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')
    """#################################################################################"""    

    ##################################   NVLM EXPORT     ##################################
    with open(f'{output_dir}/{pts_filename}.pts', 'w') as file:
        file.write('######## Panel parameters ########\n')
        file.write(f'type={surf_type}\n')
        file.write(f'n_blades={n_blades}\n\n')
        file.write(f'rotation_center={rotation_center}\n')
        file.write(f'rotation_axis={rotation_axis}\n')
        
        file.write('n_span_all= {}\n'.format(len(ALL_sections)))
        file.write('n_points={}\n'.format(len(ALL_sections[0][:,0])))
        file.write('######## End of parameters ########\n')
        
        for i in range(len(ALL_sections)):
            np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

    ########################################################################################        
    ##### Write coordinates to txt to check with the paraview ######
    with open(f'{output_dir}/Blade_points_check.txt', 'w') as file:
        file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
        for i in range(len(ALL_sections)):
            np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

    # ## Export all cross sections (Optional) ##
    # cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")

else:
    for l in range(len(collective_pitch_increment)):
        for k in range(len(twist_local_root)):
            for j in range(len(taper_local_root)):

                pts_filename_new = f'{pts_filename}_pitch_{collective_pitch_increment[l]}_twist_{twist_local_root[k]}-{twist_local_tip[k]}_taper_{taper_local_root[j]}-{taper_local_tip[j]}'

                DUST_dir = f'{working_dir}/DUST_output_pitch_{collective_pitch_increment[l]}_twist_{twist_local_root[k]}-{twist_local_tip[k]}_taper_{taper_local_root[j]}-{taper_local_tip[j]}'

                twist_local = np.linspace(twist_local_root[k],twist_local_tip[k], len(z_planes))
                taper_local = np.linspace(taper_local_root[j],taper_local_tip[j], len(z_planes))
              
                ###### GENERATE POINTS ######
                ALL_sections, ALL_sections_DUST, all_sections_compound = points.get_coords(stp_file, N_chord, dist_airfoil, z_planes, close_TE, collective_pitch_increment[l], twist_local, taper_local, find_LE)
                    
                ###### GENERATE MESH AND EXPORT #######
                mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, n_blades, close_TE)

                import meshio
                meshio.write(f'{output_dir}/{pts_filename_new}.vtk', mesh, file_format='vtk')

                """####################  DUST EXPORT (for Basic Mesh)   ###########################"""
                if os.path.isdir(DUST_dir) is False:
                    os.mkdir(f'{DUST_dir}')

                meshio.write(f'{DUST_dir}/mesh_DUST.vtk', mesh_DUST, file_format='vtk')

                with open(f'{DUST_dir}/rr.dat', 'w') as file:
                    np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

                ## Export the DUST connectivity by switching from python index to dust index
                connectivity_DUST = connectivity_DUST + np.ones(connectivity_DUST.shape)
                with open(f'{DUST_dir}/ee.dat', 'w') as file:
                    np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')
                """#################################################################################"""    

                ##################################   NVLM EXPORT     ##################################
                with open(f'{output_dir}/{pts_filename}.pts', 'w') as file:
                    file.write('######## Panel parameters ########\n')
                    file.write(f'type={surf_type}\n')
                    file.write(f'n_blades={n_blades}\n\n')
                    file.write(f'rotation_center={rotation_center}\n')
                    file.write(f'rotation_axis={rotation_axis}\n')
                    
                    file.write('n_span_all= {}\n'.format(len(ALL_sections)))
                    file.write('n_points={}\n'.format(len(ALL_sections[0][:,0])))
                    file.write('######## End of parameters ########\n')
                    
                    for i in range(len(ALL_sections)):
                        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

                ########################################################################################        
                ##### Write coordinates to txt to check with the paraview ######
                with open(f'{output_dir}/Blade_points_check.txt', 'w') as file:
                    file.write('Node_ID\tX(mm)\tY(mm)\tZ(mm)\n')
                    for i in range(len(ALL_sections)):
                        np.savetxt(file, ALL_sections[i], delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

                # ## Export all cross sections (Optional) ##
                # cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")

###############################################################################################################################
end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds")