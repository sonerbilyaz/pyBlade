from timeit import default_timer as timer
start_time = timer()  # Start the timer

from pyBlade import run
from pyBlade.in_out.config_reader import read
from pyBlade.Distributions.distributions import spanwise_disribution as spanwise_planes

###   First, read the cfg file    ###
config = read('VX4_Front_prop.cfg')

if config["GENERATE_SURFACE"]["generate_blade"] not in ['yes', True, 'Yes']:
    run.run_mesh(config)

elif config["GENERATE_SURFACE"]["generate_blade"] in ['yes', True, 'Yes']:

    run.run_blade(config)

end_time = timer()  # End the timer
print(f"Code executed in: {end_time - start_time:.6f} seconds") 