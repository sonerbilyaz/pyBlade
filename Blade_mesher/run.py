from timeit import default_timer as timer
start_time = timer()  # Start the timer

import meshio
import numpy as np, os, sys 

# Get the absolute path to the directory and add it to the sys path for relative imports
package_dir = os.path.abspath(os.path.dirname(__file__))
if package_dir not in sys.path:
    sys.path.append(package_dir)

from Distributions.distributions import spanwise_disribution as spanwise_planes
from Point_coordinates import point_generation as points
from Mesher.mesh import generate_mesh

config = {"input_dir": None, "file": None, "output_dir": None, 
          "find_LE": None, "remove_TE":None , "close_TE":None , "CCW": None,
          "N_chord": None, "N_span": None, "z_min": None, "z_max": None, 
          "dist_section": None, "dist_spanwise": None, "r_R": None,
          "collective_pitch": None, "twist_root": None, "twist_tip": None, "chord_scale": None
          }

""" ################################ INPUTS ####################################### """
#-------------------     File Properties     ------------------------
# File paths and the STEP file
config["input_dir"] = '../Runs/VX4_Aft_Rotor'
config["file"] = 'VX4_Front_Blade_single.stp'
case_prefix = 'VX4-front_prop'

# config["input_dir"] = '../Runs/VX4_Front_Prop'
# config["file"] = 'Aft_1Blade_Mesh_Orient.stp'
# case_prefix = 'VX4-aft_prop'

config["output_dir"] = f'{config["input_dir"]}/output'


## Are upper and lower surfaces seperated ???
#   If yes ==> LE is determined (find_LE = False) 
#   If NO  ==> LE is NOT determined (find_LE = True) 
config["find_LE"] = False

## Rotor Rotation direction (CCW)
#   If True (CCW)   ==> Leave x-coordinate as it is
#   If False (CW)   ==> Reverse x-coordinate direction
config["CCW"] = False

#-------------------     Panel Discretization     -------------------
## Chordwise Distribution ##
config["N_chord"] = 23                          ## Upper and Lower surf separately!!
config["dist_section"] = 'cosine_LE'

## Spanwise Distribution ##
config["N_span"]= 62
config["z_min"], config["z_max"] = 26, 150      ## in mm

config["dist_spanwise"] = 'cosine_TIP'
config["r_R"] = 0.82                            ## Span location to start the cosine_TIP 

## TE Modification ##
config["remove_TE"] = True            # Should we remove TE ??
config["close_TE"] = True             # Should we close the TE gap ??

#-----------------     Initialize section change     ------------------
### Create spanwise planes ###
z_planes = spanwise_planes(config["z_min"], config["z_max"], config["N_span"], config["dist_spanwise"], config["r_R"])

config["collective_pitch"] = np.array([0])
config["chord_scale"] = np.array([1])
config["twist_root"] = np.array([0])

config["twist_tip"] = np.array([0])

#-----------------------     Section Change     -----------------------
# ## Collective Pitch Modification (deg) ##
# config["collective_pitch"] = np.linspace(0,20,2)

# ### Spanwise local Changes ###
# config["chord_scale"] = np.linspace(0.7,1.3,2)
# config["twist_root"] = [0,10]
# config["twist_tip"] = [0,0]

#-------------------------     .pts Inputs     -------------------------
n_blades = 1                # Number of blades
surf_type = 'propeller'     # Surface type (wing/propeller)
rotation_center = [0, 0, 0]
rotation_axis = [0, -1, 0]

""" ############################################################################## """



# Create output directory if it is absent #
if os.path.isdir(config["output_dir"]) is False:
    os.mkdir(config["output_dir"])


### Create flag for the output file name ###
cond_no_change = (config["collective_pitch"] == 0).all() and (config["chord_scale"] == 1).all() and (config["twist_root"] == 0).all()


def run(config, z_planes):
    #########   GENERATE POINTS  ##########
    ALL_sections, ALL_sections_DUST, all_sections_compound = points.get_coords(config, z_planes, collective_pitch=0, twist_local=np.zeros_like(z_planes), chord_scale=1)
        
    ######### GENERATE MESH AND EXPORT #######
    mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, n_blades, config["close_TE"])

    meshio.write(f'{config["output_dir"]}/{case_prefix}.vtk', mesh, file_format='vtk')

    #########    DUST EXPORT (for Basic Mesh)   #########
    DUST_dir = f'{config["output_dir"]}/DUST_output'
    if os.path.isdir(DUST_dir) is False:
        os.mkdir(f'{DUST_dir}')

    meshio.write(f'{DUST_dir}/mesh_DUST.vtk', mesh_DUST, file_format='vtk')

    with open(f'{DUST_dir}/rr.dat', 'w') as file:
        np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

    ## Export the DUST connectivity by switching from python index to dust index
    connectivity_DUST = connectivity_DUST + np.ones(connectivity_DUST.shape)
    with open(f'{DUST_dir}/ee.dat', 'w') as file:
        np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')

    ##################################   NVLM EXPORT     ##################################
    ALL_sections_array = np.concatenate(ALL_sections, axis=0)

    with open(f'{config["output_dir"]}/{case_prefix}.pts', 'w') as file:
        file.write('######## Panel parameters ########\n')
        file.write(f'type={surf_type}\n')
        file.write(f'n_blades={n_blades}\n\n')
        file.write(f'rotation_center={rotation_center}\n')
        file.write(f'rotation_axis={rotation_axis}\n')
        
        file.write('n_span_all= {}\n'.format(len(ALL_sections)))
        file.write('n_points={}\n'.format(len(ALL_sections[0][:,0])))
        file.write('######## End of parameters ########\n')
        
        ## Write the section coordinates ##                    
        np.savetxt(file, ALL_sections_array, delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

    ########################################################################################        

def run_modify_section(config, z_planes):

    for collective_pitch in config["collective_pitch"]:
        for k in range(len(config["twist_root"])):
            for chord_scale in config["chord_scale"]:

                case_name = '{}_pitch_{:.2f}_twist_{:.2f}-{:.2f}_chord_{:.2f}'.format(case_prefix, float(collective_pitch), config["twist_root"][k] ,config["twist_tip"][k], chord_scale)
                dir_case = f'{config["output_dir"]}/{case_name}'

                ## If folder does not exist, create it
                if os.path.isdir(dir_case) is False:
                    os.mkdir(dir_case)


                twist_local = np.linspace(config["twist_root"][k], config["twist_tip"][k], len(z_planes))

                ###### GENERATE POINTS ######
                ALL_sections, ALL_sections_DUST, all_sections_compound = points.get_coords(config, z_planes, collective_pitch, twist_local, chord_scale)
                    
                ###### GENERATE MESH AND EXPORT #######
                mesh, mesh_DUST, connectivity_DUST, coordinates_DUST = generate_mesh(ALL_sections, ALL_sections_DUST, n_blades, config["close_TE"])

                meshio.write(f'{dir_case}/{case_name}.vtk', mesh, file_format='vtk')

                #########    DUST EXPORT (for Basic Mesh)   #########
                DUST_dir = f'{dir_case}/DUST_output'

                if os.path.isdir(DUST_dir) is False:
                    os.mkdir(f'{DUST_dir}')

                meshio.write(f'{DUST_dir}/mesh_DUST.vtk', mesh_DUST, file_format='vtk')

                with open(f'{DUST_dir}/rr.dat', 'w') as file:
                    np.savetxt(file, coordinates_DUST, delimiter='\t', fmt=['%.8f','%.8f','%.8f'], comments='')

                ## Export the DUST connectivity by switching from python index to dust index
                connectivity_DUST = connectivity_DUST + np.ones(connectivity_DUST.shape)
                with open(f'{DUST_dir}/ee.dat', 'w') as file:
                    np.savetxt(file, connectivity_DUST, delimiter='\t', fmt=['%.0f','%.0f','%.0f','%.0f'], comments='')

                ##################################   NVLM EXPORT     ##################################
                ALL_sections_array = np.concatenate(ALL_sections, axis=0)

                with open(f'{dir_case}/{case_name}.pts', 'w') as file:
                    file.write('######## Panel parameters ########\n')
                    file.write(f'type={surf_type}\n')
                    file.write(f'n_blades={n_blades}\n\n')
                    file.write(f'rotation_center={rotation_center}\n')
                    file.write(f'rotation_axis={rotation_axis}\n')
                    
                    file.write('n_span_all= {}\n'.format(len(ALL_sections)))
                    file.write('n_points={}\n'.format(len(ALL_sections[0][:,0])))
                    file.write('######## End of parameters ########\n')

                    ## Write the section coordinates ##                    
                    np.savetxt(file, ALL_sections_array, delimiter='\t', fmt=['%.0f','%.9f','%.9f','%.9f'], comments='')

###############################################################################################################################

if cond_no_change:
    run(config, z_planes)

else:
    run_modify_section(config, z_planes)


# ## Export all cross sections (Optional) ##
# cq.exporters.export(all_sections_compound, f"{output_dir}/all_cross_sections.step")

end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds")