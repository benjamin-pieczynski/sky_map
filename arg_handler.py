#-----------------------------------------------------------------------------------------#
# MODULE: arg_handler.py
# BY: Benjamin Pieczynski DATE: 2024-04-06
# 
# PURPOSE:
#   This module is responsible for handling the arguments provided via the command line
#   interface. The module uses a dictionary to store arguments from the command line, and
#   boolean switches.
#
# FUNCTIONS:
#   check_latlon, check_obs_time, arg_handler
#
# MODIFICATION HISTORY:
#   2024-04-06: initial development
#
#-----------------------------------------------------------------------------------------#

# imports - dependencies
import os
import time

# imports - user-defined
from utils import update_progress_bar, build_time_array
from defaults import *

# helper functions
def check_latlon(args: dict, arg_dict: dict) -> dict:
    """
        Checks the latitude and longitude for user input errors.
        Return the values by modifying the argument dictionary.
        
        In:
        args: dict - command line arguments
        arg_dict: dict - checked arguments dictionary
        
        Out:
        args_dict: dict - modified argument dictionary
    """
    try:
        arg_dict['lat'] = float(args['latitude'])
    except:
        raise ValueError('Latitude value could not be converted to Type: float')
    try:
        arg_dict['lon'] = float(args['longitude'])
    except:
        raise ValueError('Longitude value could not be converted to Type: float')
    return arg_dict

def check_user_time(obs_time: str):
    """
        Checks the observation time for user input errors.
        Return the observation values by modifying the
        argument dictionary.
        
        In:
        obs_time: str - command line arguments
        
        Out:
        None
    """
    obs_elements = obs_time.split('-')
    
    # year
    yyyy = obs_time.split('-')[0]
    if len(yyyy) != 4:
        len_year = len(yyyy)
        raise ValueError(f'Year length must be 4, got {yyyy} of length: {len_year}')
    try:
        yyyy = int(yyyy)
    except:
        raise ValueError(f'Year - {yyyy} cannot be converted to Type: Int')
    
    # month
    try:
        mm = int(obs_elements[1])
    except:
        mm = obs_elements[1]
        raise ValueError(f'Month - {mm} cannot be converted to Type: Int')
    
    # day
    try:
        dd = int(obs_elements[2][0:2])
    except:
        dd = obs_elements[2]
        raise ValueError(f'Day - {dd} cannot be converted to Type: Int')
    
    # time parsing
    hms = obs_elements[-1][2:]
    
    # check for T
    if hms[0] != 'T':
        raise ValueError('Formatting Issue. Must separarte YYYY-MM-DD and HH:MM:SS with T')
    
    # check hours
    try:
        hh = int(hms[1:3])
    except:
        hh = hms[1:3]
        raise ValueError(f'Hours - {hh} cannot be converted to Type: Int')
    
    # check minutes
    try:
        mm = int(hms[4:6])
    except:
        mm = hms[4:6]
        raise ValueError(f'Minutes - {mm} cannot be converted to Type: Int')
    
    # check seconds
    try:
        ss = int(hms[8:])
    except:
        ss = hms[8:]
        raise ValueError(f'Seconds - {ss} cannot be converted to Type: Int')

    return

def check_multitime(obs_time: str, args: dict, arg_dict: dict) -> dict:
    """
        Checks end time, time step, and units. Builds time array for
        observation times.
        
        In:
        obs_time: str - starting observation time from previous check
        args: dict - filled with command line arguments
        arg_dict: dict - inout to store values
        
        Out:
        arg_dict: dict - modified values
    """
    # set starting time
    tstart = obs_time
    
    # check ending time
    tend = args['end_time'].upper()
    check_user_time(tend)
    
    # check units
    units = args['units'].lower()
    valid_units = ['days', 'hours', 'minutes', 'seconds']
    if args['units'] not in valid_units:
        raise ValueError(f'Units - {units} must be one of {valid_units}')
    else:
        None
        
    # check time step
    try:
        h = int(args['time_step']) # prevent timestep from being float
    except:
        raise ValueError('Time Step required and must be and integer.')
    
    # all checks have passed so build the time array
    time_array = build_time_array(tstart, tend, h, units)

    # check to make sure the number of times does not exceed the limit
    len_time_array = len(time_array)
    if len_time_array > max_times:
        raise ValueError(f'Number of times ({len_time_array}) exceed program limit')
    
    # store array in dictionary
    else:
        arg_dict['obs_times'] = time_array
        
    return arg_dict

def check_paths(path: str) -> bool:
    """Purpose of this function is to check if the path
       exists so that the program can properly notify the
       user of any existing issues.\n
       INPUTS:\n
       paths: str file path
       OUTPUTS:
       bool True/False"""

    if os.path.exists(path):
        return True
    else:
        raise ValueError(f'PATH: {path} does not exist')

def check_threading(args: dict, arg_dict: dict):
    """
        Checks to see if multiprocessing is disabled.
        Reduces number of threads if the number of files is low.
        
        In:
        args: dict - argument dictionary
        
        Out:
        n_threads: int - number of threads
    """
    n_threads = int(args['n_threads'])
    bMP = args['remove_threading']
    if bMP == False:
        return 1
        
    if len(arg_dict['obs_times']) < n_threads:
        n_threads = len(arg_dict['obs_times'])
    return n_threads

# primary function
def arg_handler(args: dict) -> dict:
    
    # initialize the argument dictionary
    arg_dict = {}
    
    # progress message and tracking
    message = 'Checking Arguments'
    n_comp, n_len = 0, 6

    # check lat lon arguments
    update_progress_bar(message, n_comp, n_len, optional='checking lat lon')
    arg_dict = check_latlon(args, arg_dict)
    n_comp += 1
    time.sleep(0.2)
    
    # check timing
    update_progress_bar(message, n_comp, n_len, optional='checking times')
    obs_time = args['obs_time'].upper()
    check_user_time(obs_time)
    if args['end_time'] == False:
        arg_dict['obs_times'] = [obs_time]
    else:
        update_progress_bar(message, n_comp, n_len, optional='building time array')
        arg_dict = check_multitime(obs_time, args, arg_dict)
    n_comp += 1
    time.sleep(0.2)
    
    # check output directory 
    update_progress_bar(message, n_comp, n_len, optional='checking paths')
    out_dir = args['output_directory']
    check_paths(out_dir)
    arg_dict['out_dir'] = out_dir
    n_comp += 1
    time.sleep(0.2)
    
    # check threading
    update_progress_bar(message, n_comp, n_len, optional='checking multiprocessing')
    arg_dict['n_threads'] = check_threading(args, arg_dict)
    n_comp += 1
    time.sleep(0.2)
    
    # checking for target
    update_progress_bar(message, n_comp, n_len, optional='checking target')
    arg_dict['targets'] = args['targets']
    n_comp += 1
    time.sleep(0.2)
    
    # filling in boolean statements
    update_progress_bar(message, n_comp, n_len, optional='filling bools')
    arg_dict['bPlanets'] = args['no_planets']
    arg_dict['bSun'] = args['no_sun']
    arg_dict['bMoon'] = args['no_moon']
    arg_dict['bLabels'] = args['no_labels']
    arg_dict['bBackground'] = args['no_background']
    arg_dict['bPaths'] = args['no_paths']
    arg_dict['bTwinkling'] = args['scintillation']
    arg_dict['bTarget'] = True if args['targets'] != None else False
    arg_dict['bCoords'] = args['no_coords']
    n_comp += 1
    
    update_progress_bar(message, n_comp, n_len, optional='complete', newline=True)  
    
    return arg_dict