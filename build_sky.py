#-----------------------------------------------------------------------------------------#
# MODULE: build_sky.py
# BY: Benjamin Pieczynski DATE: 2024-04-06
# 
# PURPOSE:
#   This module is responsible for building the tables with sky coordinate positions using
#   astroquery.
#
# FUNCTIONS:
#   check_latlon, check_obs_time, arg_handler
#
# MODIFICATION HISTORY:
#   2024-04-06: initial development
#
#-----------------------------------------------------------------------------------------#

# imports
import os
import numpy as np
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_body
from astropy.table import Table, Column
from astropy.time import Time
import multiprocessing as mp
import time

# imports user-defined
from utils import update_progress_bar, convert_to_altaz, assign_threads
from plot_star_map import plot_star_map
from defaults import prog_dir

def check_hemisphere(lat: float) -> tuple:
    """
        Checks which hemisphere the observers location is
        and then returns information on the celestial pole.
        
        In:
        lat: float - latitude
        
        Out:
        cp_coords: tuple - (lon, lat, marker, celestial pole, zenith_dec)
    """
        
    if lat >= 0:
        pole = 'NCP'
        marker = '^'
        zenith_dec = 90 - lat
    elif lat < 0:
        pole = 'SCP'
        marker = 'v'
        zenith_dec = -90 - lat
    else:
        raise ValueError('Could not locate celestial pole')
    return (0, lat, 'v', pole, zenith_dec)

def query_stars(zenith_dec, mag_limit = 6.5):
    """
        Query star data from Vizier.
        
        In:
        zenith_dec
    """
        
    v = Vizier(columns=['Name', 'HD', 'RAJ2000', 'DEJ2000', 'Vmag'], 
               column_filters={"Vmag":f"<{mag_limit}"})
    v.ROW_LIMIT = -1 # set row limit to unlimited
    result = v.query_region(SkyCoord(ra=0, dec=zenith_dec, 
                                           unit=(u.deg, u.deg),
                                           frame='icrs'), 
                            radius=180*u.deg, 
                            catalog='V/50')
    
    return result[0]

def positions_to_plot(location, obs_times: list, stars, commmon_names: dict, 
                      cp_coords: tuple, arg_dict: dict, n_total, 
                      out_dir: str, qcount, lock):
    """
        Interfaces with multiprocessing to get object locations and plot
        the sky maps.
        
        In:
        location: tuple - containing EarthLocation object
        obs_times: list - array of time strings to search
        stars: astropy table - table containing star information
        common_names: dict - contains star common names
        cp_coords: tuple - celestial coordinates
        n_total: int - total number of threads running
        out_dir: str - output directory
        qcount: mp object - stores counting object
        lock: mp object - protects access of variables across threads
    """
    
    # setting update parameters
    total = len(obs_times)
    update_val  = 0.01 # update every 1% progress
    n_progress  = 0
    
    # loop through each observation time
    for obs_time in obs_times:
        
        # duplicate table
        stars_T = stars.copy()
        
        # extract star coordinates
        stars_ra  = stars['RAJ2000']
        stars_dec = stars['DEJ2000']

        # convert RA and Dec to Alt-Az coordinates
        stars_alt, stars_az = convert_to_altaz(stars_ra, stars_dec, location, obs_time)
        
        # add columns to the table
        if 'AZ' in stars_T.colnames:
            stars_T.replace_column('AZ', Column(stars_az,  name='AZ'))
        else:
            stars_T.add_column(Column(stars_az,  name='AZ'))
        
        if 'ALT' in stars_T.colnames:
            stars_T.replace_column('ALT', Column(stars_alt, name='ALT'))
        else:
            stars_T.add_column(Column(stars_alt, name='ALT'))
        
        # plot star map
        plot_star_map(stars_T, cp_coords, commmon_names, obs_time, 
                      location, arg_dict, out_dir, lock)
        
        # check to see if progress bar should be updated
        n_progress += 1
        t_progress = n_progress / n_total
        if t_progress >= update_val :
            lock.acquire()
            qcount.value += n_progress
            count = qcount.value
            lock.release()
            
            # update progress bar, reset progress
            update_progress_bar('Building Sky Maps', count, n_total)
            n_progress = 0
        else:
            None
    
    # after complete update progress bar 
    lock.acquire()
    qcount.value += n_progress
    count = qcount.value
    lock.release()
    update_progress_bar('Building Sky Maps', count, n_total)

    return

def build_sky_maps(lat: float, lon: float, arg_dict: dict, obs_times: list, n_threads: int, out_dir: str):
    """
        Builds tables with sky positions for local objects.
        
        In:
        lat: float - lattitude in degrees
        lon: float - longitude in degrees
        arg_dict: dict - with boolean arguments
        obs_times: list - array of observation times
        n_threads: int - number of threads for multiprocessing
        out_dir: str - output directory
    """
    # set user location
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
    
    # check user hemisphere
    cp_coords = check_hemisphere(lat)
    
    # query star data from Vizier
    print('querying catatlog - V/50...')
    stars = query_stars(cp_coords[-1])
    
    # get star ra and dec
    stars_ra  = stars['RAJ2000']
    stars_dec = stars['DEJ2000']

    
    # read IAU common name catalog
    filename = os.path.join(prog_dir, 'prog_data', 'IAU-CSN.csv')
    iau_table = Table.read(filename, format='csv')
    
    # build IAU dictionary for stars brighter than 2.5 mag
    print('finding common star names...')
    common_names = {}
    for row in iau_table:
        try:
            bVmag = float(row['Vmag']) < 2.7
            common_names[int(row['HD'])] = row['Name'] if bVmag else None
        except:
            None
    
    # assign observation times to threads
    print('assigning threads...')
    thread_assignments = assign_threads(obs_times, n_threads)
    
    # set up multiprocessing
    manager = mp.Manager()
    qcount = manager.Value('i', 0)
    lock = manager.Lock()
    update_progress_bar('Building Sky Maps', 0, 1)
    
    # start processes
    time_start = time.time()
    processes = []
    n_total = len(obs_times)
    for thread_key in thread_assignments:
        t_obs_times = thread_assignments[thread_key]
        p = mp.Process(target=positions_to_plot, args=(location, t_obs_times, stars, 
                                                       common_names, cp_coords, arg_dict, 
                                                       n_total, out_dir, qcount, lock))
        processes.append(p)
        p.start()
    
    # join processes
    for p in processes:
        p.join()
    
    # update progress bar
    update_progress_bar('Building Sky Maps', 1, 1, 
                        optional='Complete', newline=True)
    time_stop = time.time()
    time_total = round(time_stop - time_start, 2)
    print(f'TIME: {time_total} seconds')
    
    return