#-----------------------------------------------------------------------------------------#
# MODULE: cli.py
# BY: Benjamin Pieczynski DATE: 2024-04-06
# 
# PURPOSE:
#   This module is responsible for organizing the command line interface with the program.
#   To accomplish this I use the built in module argparse to handle argument parsing.
#
# MODIFICATION HISTORY:
#   2024-04-06: initial development
#
#-----------------------------------------------------------------------------------------#

# imports - dependencies
import os
import argparse

# imports - user-defined
from defaults import *

# program description
description = f"""
                  PROGRAM: {prog_name}
                  VERSION: {version}
                  
                  BY: {programmer} - {release_date}
                  
                  This program builds a skymap given a user location and time. 
                  The map is presented as a .png file with Alt-Az coordinates. 
                  Multiple options are provided so that different objects can 
                  be tracked or options for the object appearing on the map.
                """

# set up parser
parser = argparse.ArgumentParser(prog='sky_maps', 
                                 description=description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

# help information
lat_help = '''Observer latitude in degrees'''
lon_help = '''Observer longitude in degrees'''
obs_time = '''Observation time in format: YYYY-MM-DDTHH:MM:SS 
              (example: 2024-01-30T00:00:00)'''
out_dir_help = '''Output directory for image. Default is the current 
                  working directory'''
thread_help = '''Number of threads for multiprocessing. Default is the number of cores on
                 the current machine.'''
remove_mp_help = '''This flag removes the multiprocessing functionality. Equivalent to
                    "-nt 1"'''
endtime_help = '''Provides an end time for building multiple sky maps over a time range.
                  Stepsize will result in a endtime <= this time.'''
timestep_help = '''Time step that marches towards the end time. Default unit is seconds.'''
unit_help = '''Units for the time step. Can set to "days", "hours", "minutes", 
               or "seconds" (default).'''
no_planets_help = '''Removes planets from the plot.'''
no_sun_help  = '''Removes the Sun from the plot.'''
no_moon_help = '''Removes the Moon from the plot.'''
no_labels_help = '''Removes the object labels from the plot.'''
no_background_help = '''Removes the background from the plot.'''
target_help = '''Provide a target to track. Provides both location and path. Input as "mars venus ... saturn"'''
no_paths_help = '''Removes the paths for the Sun and Moon.'''
scintillation_help = '''Adds star scintillation, for star brightness variation between images.'''

# positional arguments
parser.add_argument('latitude', help=lat_help )
parser.add_argument('longitude', help=lon_help)
parser.add_argument('obs_time',  help=obs_time)

# optional arguments
parser.add_argument('-et', '--end_time', default=False, help=endtime_help)
parser.add_argument('-hs', '--time_step', default=3, help=timestep_help)
parser.add_argument('-u', '--units', default='seconds', help=unit_help)
parser.add_argument('-od', '--output_directory', default=cwd, help=out_dir_help)
parser.add_argument('-nt', '--n_threads', default=default_threads, help=thread_help)
parser.add_argument('-rt', '--remove_threading', action='store_false', help=remove_mp_help)
parser.add_argument('-np', '--no_planets', action='store_false', help=no_planets_help)
parser.add_argument('-ns', '--no_sun',  action='store_false', help=no_sun_help )
parser.add_argument('-nm', '--no_moon', action='store_false', help=no_moon_help)
parser.add_argument('-nl', '--no_labels', action='store_false', help=no_labels_help)
parser.add_argument('-nb', '--no_background', action='store_true', help=no_background_help)
parser.add_argument('-ne', '--no_paths', action='store_false', help=no_paths_help)
parser.add_argument('-sc', '--scintillation', action='store_true', help=scintillation_help)
parser.add_argument('-t', '--targets', nargs='*', default=None, help=target_help)

# I want to add more options for object paths and plot options