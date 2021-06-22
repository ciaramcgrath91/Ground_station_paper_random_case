"""
Created by Ciara McGrath 29/03/2021

* DISCLAIMER *
  The S/W remains property of Ciara McGrath and shall not be modified,
  nor shall it be forwarded to third parties without prior written consent

"""

import numpy as np          # need to import modules to be used in every individual module that calls it
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time
import random
from support_funcs import (
                                              odeFixedStepnoman, greg2jd, gast, eq_of_motion_noman, find_sat_elevation, 
                                              total_time_in_view, lat_long_ssp_calc, find_dist 
                                              )
                                              
def divisionbyzeroiszero(n, d):
    return n / d if d else 0

def generate_random_sources(num_sources):
    sources = []
    count = 0
    while count < num_sources:     # loop over sources
        source_temp = []                       # create temporary holder
        u = random.random()
        v = random.random()
        lat_temp = np.arccos(2*v - 1) - np.pi/2           # generate random latitude (-pi to pi) (radians) (https://mathworld.wolfram.com/SpherePointPicking.html)
        lon_temp = 2 * np.pi * u                # generate random longitude (radians) (https://mathworld.wolfram.com/SpherePointPicking.html) 
        source_temp.append(np.degrees(lat_temp))      # append random latitude in degrees
        source_temp.append(np.degrees(lon_temp))    # append random longitude in degrees
        sources.append(source_temp)                             # append to main list
        count += 1
    return sources

def generate_random_satellites(num_satellites):
    # create list of lists of satellite parameters at epoch in keplerian elements:
    # [a, e, i, raan, arg p, mean anomaly] (all in radians)
    satellites = []
    count = 0
    while count < num_satellites:     # loop over satellites
        satellite_temp = []                       # create temporary holder
        # calculate parameters
        alt = random.uniform(400, 1000)
        e = 0                               # circular orbits for now
        incl = random.uniform(0, 180)
        raan = random.uniform(0, 360)
        aop = random.uniform(0, 360)
        mean_anom = random.uniform(0, 360)
        # create list of variables in radians
        satellite_temp.append(alt*1000 + Re)      # append random latitude in degrees
        satellite_temp.append(e)    # append random longitude in degrees
        satellite_temp.append(np.radians(incl))    # append random longitude in degrees
        satellite_temp.append(np.radians(raan))    # append random longitude in degrees
        satellite_temp.append(np.radians(aop))    # append random longitude in degrees
        satellite_temp.append(np.radians(mean_anom))    # append random longitude in degrees
        # add satellite to list of satellites
        satellites.append(satellite_temp)                             # append to main list
        count += 1                                  # increment
    return satellites



start_time = time.time()

case_name = "random_test"

## ---------------------------
## Input Values
## ---------------------------

## Constants
mu =  3.986004418 * 10**14    # standard gravitational parameter, m^3/s^2
Re = 6371000               # mean Earth radius, m
J2 = 1082.7 * (10**-6)        # coefficient of the Earth's gravitational zonal harmonic of the 2nd degree, -
flattening = 0.00335281    # flattening value
rot_rate = 0.0000729212   # rotation rate of central body

## Parameters
max_target_elevation = 15 # degs
max_gstation_elevation = 15 # degs
max_ISL_range = 0 # km

## Spacecraft set up
"""
To generate randomly
"""
num_sats = 100
#sats_epoch = generate_random_satellites(num_sats)
# save satellites so can repeat analysis
#np.save("Data\\satellites_" + case_name, sats_epoch)
"""
load previously saved satellites so can replicate results
"""
sats_epoch = np.load("Data\\satellites_random_test.npy")
num_sats = len(sats_epoch)                              # number of sources  
sats_names = ["Sat " + str(i + 1) for i in range(0,num_sats)]

## Sources set up
"""
To generate randomly
"""
# Want to generate randomly
# set up number to generate
#num_sources = 300
#sources = generate_random_sources(num_sources)
# save sources so can repeat analysis
#np.save("Data\\sources_" + case_name, sources)
"""
To use previously generated sources
"""
# load previously saved sources so can replicate results
sources = np.load("Data\\sources_" + case_name + ".npy")
num_sources = len(sources)                              # number of sources   
source_names = ["Source " + str(i + 1) for i in range(0,num_sources)]


## ground stations set up: changed to test specific set
gstations = np.load("Data\\GS_locations.npy").tolist()  # load in ground stations locations in latitude and longitude
num_gstations = len(gstations)                                  # number of ground stations   
gstations_names = ["Gstation " + str(i + 1) for i in range(0,num_gstations)]

## creating lists of headings for assigning to adj matrix
headings = [source_names + sats_names + gstations_names]


## date set up
"""
date_d = 1              # start day
date_m = 1              # start month
date_y = 2020           # start year
jdate_start = greg2jd(date_m, date_d, date_y)     # convert start date to Julian
"""
jdate_start = 2458701.72286614

## analysis set up 
tmin = 0                # analysis start time
num_days = 100            # number of days to run analysis
tmax = num_days * 24 * 60 * 60     # analysis end time, secs
tstep = 30              # time step, secs

## ---------------------------
## Execution
## ---------------------------

# set start of satellite count
sat_index = 0

# create list to hold orbit parameters and lat and long values for each satellite
sat_parameters = []
data_buffer_allsats = []            # holder for data on satellite
delivered_data_allsats = []         # holder for data delievered to the ground

# loop over satellites:
for satellite in sats_epoch[sat_index:]:
    y0 = [satellite[0], satellite[3], satellite[4] + satellite[5]]     # assign sma, raan and arg of latitude
    incl = satellite[2]                                 # assign inclination

    # integrate equations. sol has states of each value at all times in t vector 
    sol, t = odeFixedStepnoman(eq_of_motion_noman, y0, tmin, tmax, tstep, mu, Re, J2, incl)

    # create timeline for satellite 
    timeline = [[i] for i in t]

    # assign values
    a_full = sol[:, 0]
    raan_full = sol[:, 1]
    aol_full = sol[:, 2]

    # initialise value holders
    lat_ssp = []
    lon_ssp = []

    # loop through solutions and calculate latitude and longtude of ssp at each time step
    for i in range(0,len(sol)):
            sat_ssp_temp = lat_long_ssp_calc(sol[i], t[i], Re, incl, flattening, rot_rate, jdate_start)
            lat_ssp.append(sat_ssp_temp[0])
            lon_ssp.append(sat_ssp_temp[1])

    # set start of sources count
    source_index = 0   
    
    # loop over targets
    for source in sources:
        # set up target
        source_lat_deg = source[0]         # latitude of target, deg
        source_lon_deg = source[1]           # longitude of target, deg
        source_lat = source_lat_deg * np.pi / 180      # latitude of target, deg
        source_lon = source_lon_deg * np.pi / 180        # longitude of target, deg

        # initialise value holders
        dist = []
        elevation = []

        # loop through solutions and calculate elevation at each time step
        for i in range(0,len(lat_ssp)):
            dist_temp = find_dist(lat_ssp[i], lon_ssp[i], source_lat, source_lon, Re)
            elev_temp = find_sat_elevation(sol[i][0], dist_temp, Re)
            dist.append(dist_temp)
            elevation.append(elev_temp) 

        # convert elevation to degrees
        elevation_deg = [i * 180 / np.pi for i in elevation]

        # want instances with elevation >x deg
        el_deg_in_view = [i for i in elevation_deg if i >= max_target_elevation]
        index_in_view = [i for i, v in enumerate(elevation_deg) if v >= max_target_elevation]       # returns index i if value v is greater than or equal to max elevation: https://stackoverflow.com/questions/13717463/find-the-indices-of-elements-greater-than-x
        timesteps_in_view = [t[i] for i in index_in_view]

        # indicate in timeline when source is in view and which source
        for i in index_in_view: 
            timeline[i].append(source_names[source_index])
        
         # then need to get total time in view for each window
        if timesteps_in_view != []:     # but only if there actually are any sightings
            t_in_view = total_time_in_view(timesteps_in_view, tstep)

        source_index += 1                                                 # increment ground station count for next loop


    # set start of ground stations count
    gstation_index = 0   
    
    # loop over ground stations
    for station in gstations:
        # set up ground station
        tar_lat_deg = station[0]         # latitude of target, deg
        tar_lon_deg = station[1]           # longitude of target, deg
        tar_lat = tar_lat_deg * np.pi / 180      # latitude of target, deg
        tar_lon = tar_lon_deg * np.pi / 180        # longitude of target, deg

        # initialise value holders
        dist = []
        elevation = []

        # loop through solutions and calculate elevation at each time step
        for i in range(0,len(lat_ssp)):
            dist_temp = find_dist(lat_ssp[i], lon_ssp[i], tar_lat, tar_lon, Re)
            elev_temp = find_sat_elevation(sol[i][0], dist_temp, Re)
            dist.append(dist_temp)
            elevation.append(elev_temp) 

        # convert elevation to degrees
        elevation_deg = [i * 180 / np.pi for i in elevation]

        # want instances with elevation > x deg
        el_deg_in_view = [i for i in elevation_deg if i >= max_gstation_elevation]
        index_in_view = [i for i, v in enumerate(elevation_deg) if v >= max_gstation_elevation]       # returns index i if value v is greater than or equal to 10: https://stackoverflow.com/questions/13717463/find-the-indices-of-elements-greater-than-x
        timesteps_in_view = [t[i] for i in index_in_view]

        # indicate in timeline when sink is in view and which gstation
        for i in index_in_view: 
            timeline[i].append(gstations_names[gstation_index])
        

        # then need to get total time in view for each window
        if timesteps_in_view != []:     # but only if there actually are any sightings
            t_in_view = total_time_in_view(timesteps_in_view, tstep)

        gstation_index += 1                                                 # increment ground station count for next loop

    # save timelines
    np.save("Timeline_files\\timeline_" + "sat_" + str(sat_index) + '_' + case_name, timeline)

    sat_index +=1                           # increment ground station count for next loop


runtime = (time.time() - start_time)

print("--- runtime: %s seconds ---" % (runtime))


# save run times in text file. Appends so keep all data
file1 = open("runtimes.txt","a") #append mode 
file1.write("Runtime: " + str(runtime) + " secs for " + str(num_days) + " days, " + str(num_sources) + " sources, " + str(num_sats) + " satellites, and " + str(num_gstations) + " ground stations. \n") 
file1.close() 


