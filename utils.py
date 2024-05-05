#-----------------------------------------------------------------------------------------#
# MODULE: arg_handler.py
# BY: Benjamin Pieczynski DATE: 2024-04-06
# 
# PURPOSE:
#   This module is responsible for holding utility functions that may server the main
#   program.
#
# FUNCTIONS:
#   update_progress_bar
#
# MODIFICATION HISTORY:
#   2024-04-06: initial development
#
#-----------------------------------------------------------------------------------------#

# imports - dependencies
import os
from datetime import datetime, timedelta, timezone
from astropy.coordinates import SkyCoord, AltAz, get_sun, get_body
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.time import Time
import numpy as np



def update_progress_bar(message: str, n_comp: int, n_len: int, optional: str ='',
                        reg: bool = True, newline: bool = False) -> None:
    """
        Update the user with a print statement that displays a progress bar.
        
        In:
        message: str - program message
        n_comp: int - number completed
        n_len: int length of the array that must be processed.
        optionl: str - additional progress message
        reg: bool - whether to return register to start of line.
        newline: bool - if True, after this update go to the next line.
        
        Out:
        None
    """

    # check progress and build progress statement
    pid      = os.getpid()
    per_comp = int(100*n_comp / n_len)
    arrows   = '>' * (per_comp // 5)
    arrows   = arrows.ljust(20)
    per_comp = str(per_comp).rjust(3)
    prog_s   = f'({pid}) {message}:[{arrows}] progress - {per_comp}% {optional}'
    
    # if register return location
    if reg == True:
        print('\r' + ' ' * 100, end='')  # Clear previous line
        print('\r' + prog_s, end='')
    else:
        print(prog_s)
    
    if newline == True:
        print('\n')
    return



def build_time_array(tstart: str, tend: str, h: int, units: str) -> list:
    """
        Builds the time array from tstart to tend using h (step size)
        
        In:
        tstart: str - YYYY-MM-DDTHH:MM:SS
        tend: str - YYYY-MM-DDTHH:MM:SS
        h: int - step size
        units: str - unit name ['days', 'hours', 'minutes', 'seconds']

        Out:
        time_array: list - array of times including the start time
    """
    # parse start and end times in UTC
    start_time = datetime.strptime(tstart, "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc)
    end_time = datetime.strptime(tend, "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc)
    
    # initialize the result list with the start time
    time_array = [tstart]
    current_time = start_time
    
    # select delta
    if units == 'days':
        delta = timedelta(days=h)
    elif units == 'hours':
        delta = timedelta(hours=h)
    elif units == 'minutes':
        delta = timedelta(minutes=h)
    elif units == 'seconds':
        delta = timedelta(seconds=h)
        
    # increment time by 'h' units until reaching or exceeding the end time
    while True:
        current_time += delta
        if current_time <= end_time:
            time_array.append(current_time.strftime("%Y-%m-%dT%H:%M:%S"))
        else:
            break
    
    return time_array



def convert_to_altaz(ra, dec, location, obstime: str) -> tuple:
    """
        Converts ra and dec to alt az coordinates for a given location
        and observation time.
        
        In:
        ra: right ascension
        dec: declination
        location: lat lon location from EarthLocation
        obstime: str - format YYYY-MM-DDTHH:MM:SS
        
        Out:
        alt_deg - altitude in degrees
        az_deg - azimuth in degrees
    """
    # Convert RA and Dec from masked columns to arrays
    try:
        ra_str = ra.filled().data.tolist()
        dec_str = dec.filled().data.tolist()

        ra_deg = []
        dec_deg = []

        # Convert RA and Dec strings to degrees
        for ra_val, dec_val in zip(ra_str, dec_str):
            ra_components = ra_val.split()
            ra_hour, ra_min, ra_sec = map(float, ra_components)

            dec_components = dec_val.split()
            dec_deg_val, dec_min, dec_sec = map(float, dec_components)

            ra_deg.append((ra_hour + ra_min / 60 + ra_sec / 3600) * 15)  # Convert RA to degrees
            dec_deg.append(dec_deg_val + dec_min / 60 + dec_sec / 3600)  # Convert Dec to degrees
    except:
        ra_deg = ra
        dec_deg = dec
        
    # Convert time to the required format
    try:
        obstime = obstime.split('T')[0] + ' ' + obstime.split('T')[1]
    except:
        None # already in the correct format
    
    # Create SkyCoord object
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
   
  # Convert to Alt-Az coordinates
    altaz = coord.transform_to(AltAz(obstime=obstime, location=location))
    
    # Extract Alt-Az coordinates
    alt_deg = altaz.alt.degree
    az_deg = altaz.az.degree
    
    return alt_deg, az_deg



def assign_threads(arr: list, n_threads: int) -> dict:
    """
        Takes in a list and divides it amond threads
        
        In:
        array: list - input array
        n_threads: int - the number of threads
        
        Out:
        thread_assignments: dict - list of assignments per thread
    """
    
    # if multiprocessing is not used
    if n_threads == 1:
        return {1: arr}
    
    # thread initialization
    thread_assignments = {}
    chunk_sizes = [] # size of each thread
    
    # chunk size parameters
    chunk_size = len(arr) // n_threads  # Calculate the chunk size
    remainder = len(arr) % n_threads
    thread_num = 1
    
    # assign contents to each thread
    start = 0
    for i in range(n_threads):
        # calculate the end index of the chunk
        end = start + chunk_size + (1 if i < remainder else 0)
        thread_assignments[thread_num] = arr[start:end]
        start = end
        thread_num += 1
    
    return thread_assignments



def time_proc(obstime: str):
    """
        Builds time arrays for use in plotting object locations over time

        In:
        obstime: str
        
        Out:
        t_current: curent time as a Time object
        times_str: iso version of times
    """
    # Calculate half day times
    half_day = 0.5*u.day
    t_current = Time(obstime)
    t_start = t_current - half_day
    t_end = t_current + half_day
    
    # Generate an array of times within the specified range
    times = np.linspace(t_start.jd, t_end.jd, 1000)
    t_range = Time(times, format='jd', scale='utc')
    
    # Also create and array as a string version
    times_str = t_range.iso
    
    return t_current, times_str



def solar_system_target_pos(obstime:str, location: tuple, target: str, path=False) -> tuple:
    """
        Grabs the current position for a target and returns the translated coordinates.
        Can return path or it can only return the current coordinated.
        
        In:
        obstime: str - current time
        location: tuple - observer location in lat, lon
        target: str - name of the target object
        
        Out:
        body_coords: tuple - az, 90-alt converted coordinate positions
    """
    
    t_current = Time(obstime)    
    
    # Get the objects current position
    current_pos = get_body(target, t_current)
    body_ra  = current_pos.ra.deg
    body_dec = current_pos.dec.deg
    
    # Convert current position to Alt-Az
    body_alt, body_az = convert_to_altaz(body_ra, body_dec, location, obstime)
    body_coords = (body_az, 90 - body_alt)
    
    # Convert to Alt-Az for all the times (path tracing)
    if path:
        _, times_str = time_proc(obstime)
        path_coords = ([], []) # AZ, ALT
        for _, tobs in enumerate(times_str):
            path_alt, path_az = convert_to_altaz(body_ra, body_dec, location, tobs)
            path_coords[0].append(path_az)
            path_coords[1].append(90 - path_alt)

        return body_coords, path_coords

    else:
        return body_coords



# Function to get the ecliptic path
def ecliptic_pos(obstime: str, location: tuple):
    
    # Get the time arrays for the past 24 hours
    t_current, times_str = time_proc(obstime)
    
    # Get the sun's position at the current_time
    current_pos = get_sun(t_current)
    sun_ra_deg = current_pos.ra.deg
    sun_dec_deg = current_pos.dec.deg
    
    # Convert to Alt-Az for current sun position
    sun_alt, sun_az = convert_to_altaz(sun_ra_deg, sun_dec_deg, location, obstime)
    sun_coords = (sun_az, 90 - sun_alt)
 
    # Convert to Alt-Az for all the times
    ecliptic_coords = ([], []) # AZ, ALT
    for _, tobs in enumerate(times_str):
        ecliptic_alt, ecliptic_az = convert_to_altaz(sun_ra_deg, sun_dec_deg, location,tobs)
        ecliptic_coords[0].append(ecliptic_az)
        ecliptic_coords[1].append(90 - ecliptic_alt)
    
    return sun_coords, ecliptic_coords

def target_pos(obstime: str, location: tuple, target: str, path: bool = False) -> tuple:
    """
    Grabs the target position by checking the solar system ephemeris or by using SIMBAD.
    First it checks the solar system objects, next it checks SIMBAd to get RA and DEC.
    Then it uses the observation time to calculate the positions.
    """
    
    try:
        body_coords, body_path = solar_system_target_pos(obstime, location, target, path=True)
    except:
        try:
            t_current = Time(obstime)
            simbad = Simbad()
            simbad.add_votable_fields('ra(d)', 'dec(d)')
            result = simbad.query_object(target)
            ra_deg  = result['RA_d' ][0]  # RA in decimal degrees
            dec_deg = result['DEC_d'][0]  # DEC in decimal degrees
            
            # Convert current position to Alt-Az
            body_alt, body_az = convert_to_altaz(ra_deg, dec_deg, location, obstime)
            body_coords = (body_az, 90 - body_alt)

            # Convert to Alt-Az for all the times (path tracing)
            if path:
                _, times_str = time_proc(obstime)
                path_coords = ([], []) # AZ, ALT
                for _, tobs in enumerate(times_str):
                    path_alt, path_az = convert_to_altaz(ra_deg, dec_deg, location, tobs)
                    path_coords[0].append(path_az)
                    path_coords[1].append(90 - path_alt)
            else:
                return body_coords
        except:
            raise ValueError(f'Could not find the target using Simbad or get_body: ({target})')
            
    return body_coords, path_coords