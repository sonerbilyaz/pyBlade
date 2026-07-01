from timeit import default_timer as timer
start_time = timer()  # Start the timer


from pyBlade import run
from pyBlade.in_out.config_reader import read
from pyBlade.Distributions.distributions import spanwise_disribution as spanwise_planes

###   First, read the cfg file    ###
config = read('VX4_Front_prop.cfg')

### Create spanwise planes ###
z_planes = spanwise_planes(config["PANEL"]["z_min"], config["PANEL"]["z_max"], config["PANEL"]["N_span"], config["PANEL"]["dist_span"], config["PANEL"]["r_R"])

if config["MODIFY"]["modify_planform"] not in ['yes', True, 'Yes']:
    run.run(config, z_planes)

elif config["MODIFY"]["modify_planform"] in ['yes', True, 'Yes']:
    run.run_modify_planform(config, z_planes)

end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds")