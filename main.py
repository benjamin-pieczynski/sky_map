#-----------------------------------------------------------------------------------------#
# PROGRAM: Sky Maps
# BY: Benjamin Pieczynski DATE: 2024-04-06
# 
# DESCRIPTION:
#   This program builds a skymap given a user location and time. The map is presented as
#   a .png file with Alt-Az coordinates. Multiple options are provided so that different
#   objects can be tracked or options for the object appearing on the map.
#
# MODIFICATIONS:
#   2024-04-06: initial development
#
#-----------------------------------------------------------------------------------------#

# imports - dependencies
import numpy as np
import matplotlib.pyplot as plt
import random

# imports - user-defined
from defaults import *
from cli import parser
from arg_handler import arg_handler
from build_sky import build_sky_maps

# main program
def main():
    # print program
    print(f'\nPROGRAM: {full_name}')
    print(f'BY: {programmer}\n')
    print('\nstarting program...\n')
    # use parser to get arguments
    args = parser.parse_args().__dict__
    
    # pass arguments to argument handler for checking
    args = arg_handler(args)
    
    # store variables as needed
    obs_times = args['obs_times']
    lat = args['lat']
    lon = args['lon']
    n_threads = args['n_threads']
    out_dir = args['out_dir']
    
    print('N_THREADS =', n_threads)
    
    # build the skymaps
    build_sky_maps(lat, lon, args, obs_times, n_threads, out_dir)
    return
    
# run the program
if __name__ == '__main__':
    main()