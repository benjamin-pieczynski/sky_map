#-----------------------------------------------------------------------------------------#
# MODULE: defaults.py
# BY: Benjamin Pieczynski DATE: 2024-04-06
# 
# PURPOSE:
#   This module is responsible for setting the default settings and variables for the
#   program.
#
# MODIFICATION HISTORY:
#   2024-04-06: initial development
#
#-----------------------------------------------------------------------------------------#

# imports - dependencies
import os
import multiprocessing as mp

# meta data
prog_name = 'Sky Maps'
version = 'v0.0.0'
programmer = 'Benjamin Pieczynski'
release_date = '2024-04-06'
full_name = f'{prog_name} - {version} - {release_date}'

# globals
cwd = os.getcwd() # grabs the current working directory
default_threads = mp.cpu_count() # grabs the number of available threads
max_times = 3000 # maximum number of files that will be generated
planet_dict = {
               'Mercury': ['grey', 35], 
               'Venus': ['khaki', 70], 
               'Mars': ['orangered', 45], 
               'Jupiter': ['orange', 65], 
               'Saturn': ['goldenrod', 60], 
               'Uranus': ['darkturquoise', 30], 
               'Neptune': ['dodgerblue', 28], 
               }
