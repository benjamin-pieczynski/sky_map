#-----------------------------------------------------------------------------------------#
# MODULE: plot_star_map.py
# BY: Benjamin Pieczynski DATE: 2024-04-06
# 
# PURPOSE:
#   This module is responsible for plotting the star map.
#
# FUNCTIONS:
#   plot_star_map
#
# MODIFICATION HISTORY:
#   2024-04-06: initial development
#   2024-04-21: continued development
#
#-----------------------------------------------------------------------------------------#

# imports - dependencies
import os
import numpy as np
import random
import matplotlib.pyplot as plt
from datetime import datetime
import time

# imports - user
from utils import ecliptic_pos, solar_system_target_pos, target_pos
from defaults import planet_dict

def plot_star_map(stars, cp_coords: tuple, common_names: tuple, 
                  obs_time: str, location: tuple, arg_dict: dict,
                  out_dir: str, lock):
    
    # Grab alt-az positions and brightness
    stars_az  = stars['AZ']
    stars_alt = stars['ALT']
    stars_name = stars['Name']
    hd_ids = stars['HD']
    stars_Vmag = stars['Vmag']
    
    # Convert star brightness to sizes (adjust the scaling factor as needed)
    star_sizes = 10*(stars_Vmag*-1 + 7 + 0.01)
    
    # scintillation using alpha
    if arg_dict['bTwinkling']:
        star_alphas = [random.uniform(0.35, 1) for _ in range(len(stars_alt))]
    else:
        star_alphas = [1 for _ in range(len(stars_alt))]
    
    # plot figure
    fig = plt.figure(figsize=(60, 60))
    ax = fig.add_subplot(111, polar=True)
    
    # plot stars
    ax.scatter(np.radians(stars_az), 90 - stars_alt, color='white', s=star_sizes, 
               alpha=star_alphas, marker='o', zorder=6)
    
    # name offset
    offset_rad = 0.42 # deg
    offset_az = np.radians(0.2) # radians
    
    # search for common names
    if arg_dict['bLabels']:
        for _, (az, alt, hd_id) in enumerate(zip(stars_az, stars_alt, hd_ids)):
            try:
                name = common_names[hd_id]
                ax.annotate(name, xy=(np.radians(az) + offset_az, 90 - alt + offset_rad), 
                            color='skyblue', zorder=7)
            except:
                None
            
    # grab the sun and ecliptic coordinates for the current time.
    if arg_dict['bSun']:
        sun_coords, ecliptic_coords = ecliptic_pos(obs_time, location)
        
        # check if ecliptic option selected
        if arg_dict['bPaths']:
            
            # plot ecliptic
            ax.plot(np.radians(ecliptic_coords[0]), ecliptic_coords[1], lw=5, color='skyblue', 
                    ls='dashdot', alpha=0.5, zorder=2)
    
        # plot the Sun
        ax.scatter(np.radians(sun_coords[0]), sun_coords[1], s=300, color="gold", 
                   alpha=0.7, zorder=7)
        
        # Sun label
        if arg_dict['bLabels']:
            ax.annotate('Sol', xy=(np.radians(sun_coords[0]) + offset_az, sun_coords[1] + offset_rad), color='skyblue', zorder=7)
    
    # Plot the Moon
    if arg_dict['bMoon']:

        # grab moon coordinates
        if arg_dict['bPaths']:
            moon_coords, moon_path = solar_system_target_pos(obs_time, location, 'moon', path=True)
        else:
            moon_coords = solar_system_target_pos(obs_time, location, 'moon')
        
        # if path option plot the moon path
        if arg_dict['bPaths']:
            ax.plot(np.radians(moon_path[0]), moon_path[1], lw=5, color='lightsteelblue',
                    ls='dotted', alpha=0.5, zorder=2)
        ax.scatter(np.radians(moon_coords[0]), moon_coords[1], s=290, color='whitesmoke',
                   alpha=0.7, zorder=8)
        
        # moon label
        if arg_dict['bLabels']:
            ax.annotate('Luna', xy=(np.radians(moon_coords[0]) + offset_az, moon_coords[1] + offset_rad), color='skyblue', zorder=8)

    # Plot planets
    if arg_dict['bPlanets']:
        for planet in planet_dict:
            planet_coords = solar_system_target_pos(obs_time, location, planet.lower())
            ax.scatter(np.radians(planet_coords[0]), planet_coords[1], 
                       s=planet_dict[planet][1], color=planet_dict[planet][0])
            if arg_dict['bLabels']:
                ax.annotate(planet, xy=(np.radians(planet_coords[0]) + offset_az, planet_coords[1] + offset_rad), color='skyblue', 
                            alpha=0.7, zorder=8)
                
    # Plot targets
    if arg_dict['bTarget']:
        for target in arg_dict['targets']:
            target_coords, target_path = target_pos(obs_time, location, str(target).lower(), path=True)
            ax.scatter(np.radians(target_coords[0]), target_coords[1], marker='o', facecolors='none',
                       edgecolors='r', s=80, zorder=15)
            ax.scatter(np.radians(target_coords[0]), target_coords[1], marker='x', color='r', s=60, zorder=15)
            if arg_dict['bPaths']:
                ax.plot(np.radians(target_path[0]), target_path[1], lw=4, color='red', ls=':', alpha=0.6)
            if arg_dict['bLabels']:
                ax.annotate(target, xy=(np.radians(target_coords[0]) + offset_az, 
                            target_coords[1] + offset_rad), color='red', zorder=15)

    # adjust parameters
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 90)
    ax.set_yticklabels([])
    ax.grid(False)
    
    # Add semi-transparent filled polygon as a background to mimic grid lines
    # Plot AZ lines
    if arg_dict['bCoords']:
        bALTlines = False # alt lines filled
        
        for az_angle in range(0, 361, 30):
            ax.plot(np.radians([az_angle, az_angle]), [10, 90], color='steelblue', alpha=0.5, zorder=2)
            if az_angle != 360:
                ax.annotate(r'${}.00\degree$'.format(az_angle), xy=(np.radians(az_angle), 90),
                            color='steelblue', alpha=0.7, zorder=3)
                ax.annotate(r'${}.00\degree$'.format(az_angle), xy=(np.radians(az_angle), 50),
                            color='steelblue', alpha=0.7, zorder=3)
                ax.annotate(r'${}.00\degree$'.format(az_angle), xy=(np.radians(az_angle), 10),
                            color='steelblue', alpha=0.7, zorder=3)
                
            # Plot coordinate values
            for alt_angle in range(0, 99, 10):
                if bALTlines == False:
                    ax.plot(np.linspace(0, 2*np.pi, 1000), np.ones(1000)*alt_angle, color='steelblue', alpha=0.5, zorder=2)
                    
                if alt_angle not in [0, 10]:
                    ax.annotate(r'${}.00\degree$'.format(90-alt_angle), xy=(np.radians(az_angle + 15), alt_angle),
                                color='steelblue', alpha=0.7, zorder=3)
            
            # alt lines have been filled so make sure it is set to true
            bALTlines = bALTlines or True
        
        # Plot celestial pole
        if cp_coords[3] == 'NCP':
            ax.scatter(cp_coords[0], 90-cp_coords[1], marker=cp_coords[2], s=100, color='red', zorder=4)
        elif cp_coords[3] == 'SCP':
            ax.scatter(cp_coords[0], 90+cp_coords[1], marker=cp_coords[2], s=100, color='red', zorder=4)
        
        # Plot zenith
        #ax.annotate(r'${}.00\degree$'.format(90), xy=(0, 0), color='steelblue', alpha=0.7, zorder=3)
        ax.scatter(0, 0, marker='x', s=150, color='steelblue', alpha=0.5, zorder=2)
    
    # Remove the graph outline
    ax.spines['polar'].set_visible(False)

    # Remove tick marks
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    
    # set figure background
    if arg_dict['bBackground']:
        ax.set_facecolor('none')  # Transparent background
    else:
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')
    
    # formatting filename
    translation_table = str.maketrans({"-": "", ":": ""})
    file_time = obs_time.translate(translation_table)
    filename = os.path.join(out_dir, "skymap_" + file_time + ".png")
    
    # write file
    lock.acquire()
    plt.savefig(filename, transparent=arg_dict['bBackground'])
    #plt.show()
    lock.release()
    plt.close()
    
    return