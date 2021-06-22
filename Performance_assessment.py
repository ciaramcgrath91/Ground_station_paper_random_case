"""
Load in timelines and analyse to see performance for data delivery with subsets of ground stations
"""


import numpy as np
import pandas as pd
import time
import random
import csv
import pickle

def find_oldest_data(data_buffer, target_names, expiry, current_timestep):
    """
    Function to find the oldest data in the data buffer dictionary
    that has not already exceeded the expiry time.
    Returns ID of data in dictionary and age of data
    """
    # set these as 0 to start
    oldest_data_source = 0
    oldest_data_time = 0
    
    for source in target_names:                # loop over each holder in the data dictionary
        if data_buffer[source]:    # if the holder is empty (i.e. == None) then skip
            age_of_data = current_timestep - data_buffer[source]   # calculate age of data as current time - collection time
            if age_of_data > oldest_data_time and age_of_data < expiry:     # if it's older than the current reference AND not older than 60 mins
                oldest_data_time = age_of_data      # update the reference time
                oldest_data_source = source         # update the oldest data ID
            pass
        pass
    return oldest_data_source, oldest_data_time

def find_least_connected(data_vol_delivered_allsources, target_names, total_data_delivered):
    """
    Function to find the least connected source
    and the volume from that source
    """
    least_connected_source = None   # set as empty to start
    data_vol_least_connected = total_data_delivered # set starting guess as the total of all collected data
    
    for source in target_names:         # for each source
        data_from_source = data_vol_delivered_allsources[source]
        if data_from_source < data_vol_least_connected:     # if the data collected from this source is less than the current reference
            data_vol_least_connected = data_from_source     # update the reference for the least connected data volume
            least_connected_source = source                 # update the least connected source
            pass
        pass
    return least_connected_source, data_vol_least_connected

"""
USER DEFINED
"""
case_type = 'max flow' # option of 'max flow' or 'consensus'

files_name = "random_test"
case_tag = "40tar_10sat"

timeline_path = "Timeline_files\\"                       # define path to load timelines                                    
results_path = "Performance_files\\"                       # define path to save results                                    






"""
Set Up
"""
start_time = time.time()
# number of days simulated
num_days = 100  #7
tmax = num_days * 24 * 60 * 60     # analysis end time, secs
#set timestep length in seconds
time_step_length = 30
#number of timesteps to consider
num_timesteps = int(num_days * 24 * 60 * 60 / time_step_length)

# max time before deleting data (secs):
max_expiry_hours = 1
max_expiry = 60 * 60 * max_expiry_hours



"""
Ground Stations
"""
# set list of gstations
if case_type == 'consensus':
    gstation_filename = "GS_CNS_" + case_tag
elif case_type == 'max flow':
    gstation_filename = "GS_MF_" + case_tag

# load in gstation numbers as lists
with open('Data\\' + gstation_filename + '.csv', newline='') as f:
    reader = csv.reader(f)
    gstation_lists = list(reader)

gstations_names_list = []
for list_gstations in gstation_lists:
    gstations_names_list.append(["Gstation " + i for i in list_gstations])

"""
Satellites
"""
satellites_filename = "satellites_" + case_tag
# load in satellite numbers as lists
with open('Data\\' + satellites_filename + '.csv', newline='') as f:
    reader = csv.reader(f)
    satellite_lists = list(reader)

satellite_names_list = []
for list_satellite in satellite_lists:
    satellite_names_list.append(["sat_" + str(int(i) - 301) for i in list_satellite])    # because my count starts from 0 for spacecraft and Ruaridh's starts from 301...

"""
Targets
"""
target_filename = "targets_" + case_tag
# load in target numbers as lists
with open('Data\\' + target_filename + '.csv', newline='') as f:
    reader = csv.reader(f)
    target_lists = list(reader)

target_names_list = []
for list_tars in target_lists:
    target_names_list.append(["Source " + i for i in list_tars])  



"""
MAIN LOOP
"""
# loop over cases
for case_num in range(0, len(gstation_lists)):      # for each set of case parameters
    case_name = 'Case_' + str(case_num) + '_' + case_type + '_' + files_name + '_' + case_tag

    # set names of each ground station for this case
    gstations_names = gstations_names_list[case_num]                                                                          # set names using selected list
    num_gstations = len(gstations_names)

    # set names of each target for this case
    target_names = target_names_list[case_num]                                                                          # set names using selected list
    num_targets = len(target_names)

    # set names of each satellite for this case
    sat_names = satellite_names_list[case_num]                                                                          # set names using selected list
    num_sats = len(sat_names)

    # performance metrics holders for each case
    data_vol_delivered_allsats = { i : 0 for i in sat_names}  # dictionary size of all sats
    data_vol_delivered_allgstations = { i : 0 for i in gstations_names} # dictionary of all ground stations
    data_vol_delivered_allsources = { i : 0 for i in target_names}  # dictionary of all sources
    delay_allsats = { i : 0 for i in sat_names}  # dictionary of all sats
    delay_allsats_mins = { i : 0 for i in sat_names}  # dictionary of all sats
    delay_allsats_alldata = []


    # Analyse timeline for each satellite we are interested in in this case:
    for sat in sat_names:
        timeline = np.load(timeline_path + 'timeline_' + sat + '_' + files_name + '.npy', allow_pickle = True)
        
        """
        Do "data routing"
        """
        # create dictionary to hold time of data collection with one key for each target
        data_buffer = { i : None for i in target_names}   # set as None when "empty"
        delivered_data = []         # holder for data delievered to the ground

        time_step_count = 0
        downlink_count = 0          # counter that increments by 1/(number of sats in view of gstation) at each step a gstation is in view. downlink occurs once it's >= 1
        for time_step in timeline[0:num_timesteps]:
            if len(time_step) > 1:         # if has more than one entry, and therefore something is in view
                for sight in time_step[1:]:     # skip timestamp
                    if sight in target_names:    # if it's a source and is in this analysis and it hasn't just been see by the same satellite
                        data_buffer[sight] = time_step[0]                    # store collection time in dictionary overwriting previous entry
                    elif sight in gstations_names:    # if it's a sink and is in this analysis
                        downlink_count += 1   # increments by 1 at each step a gstation is in view.
                        if downlink_count >= 1:
                            # check which data to downlink
                            data_to_move, delay = find_oldest_data(data_buffer, target_names, max_expiry, time_step[0])
                            # data_to_move will be zero if the buffer is empty
                            if data_to_move != 0:                        
                                delivered_data.append({'target' : data_to_move, 'gstation' : sight, 'delay' : delay})   # add downlinked data to list as dictionary with target ID and delay time
                                delivered_data.append({'target' : data_to_move, 'gstation' : sight[2], 'delay' : delay})   # add downlinked data to list as dictionary with target ID and delay time
                                data_buffer[data_to_move] = None      # clear target holder in buffer so data won't be downlinked again
                                downlink_count -= 1            # take 1 from the downlink counter
        
            time_step_count += 1

        """
        Get performance metrics
        """
        data_vol_delivered_allsats[sat] = (len(delivered_data))       # get number of delivered packets for each satellite)

        # calculate data volume delivered for each ground station
        for name in gstations_names:                                            # for each ground station
            data_down_station = len([i for i in delivered_data if i['gstation'] == name])         # calculate the volume of data downlinked to the ground stations in question and add it to the total
            data_vol_delivered_allgstations[name] += data_down_station                  # add on data downlinked by each sat to each gstation -- pull gstation name from 'name' variables

        # calculate data volume delivered for each source
        for name in target_names:                                            # for each source
            data_down_tar = len([i for i in delivered_data if i['target'] == name])     # calculate the volume of data downlinked that came from the source in question and add it to the total      
            data_vol_delivered_allsources[name] += data_down_tar                  # add on data delivered by each sat from each source -- pull name from 'name' variable

        # calculate average delay for each satellite
        if delivered_data == []:    # if no delivered data
            delay_allsats[sat] = None
            delay_allsats_mins[sat] = None
        else:
            delay_allsats[sat] = sum([i['delay'] for i in delivered_data])/len(delivered_data)              # get average of all delays
            delay_allsats_mins[sat] = (sum([i['delay'] for i in delivered_data])/len(delivered_data))/60              # get average of all delays

    # total values
    data_vol_total = sum([data_vol_delivered_allsats[sat] for sat in sat_names])  # total data delivered by all satellites combined
    # total delay is average of all delays, for all satellites that have a delay value
    delay_total = sum([delay_allsats[sat] for sat in sat_names if delay_allsats[sat]])/len([delay_allsats[sat] for sat in sat_names if delay_allsats[sat]])
    delay_average_mins = delay_total/60


    # find worst connected source and volume of data from least connected source
    least_connected_source, data_vol_least_connected = find_least_connected(data_vol_delivered_allsources, target_names, data_vol_total)
    
    # Save results
    file1 = open(results_path + 'Performance_' + case_name + '_' + str(num_days) +'_days.txt',"w")
    file1.write("\nData volume delivered: " + str(data_vol_total) + "\nData volume delivered by satellite: " + str(data_vol_delivered_allsats) + "\nData volume per ground station: " + str(data_vol_delivered_allgstations) + "\nData volume downlinked per source: " + str(data_vol_delivered_allsources) + "\nMinimum data volume from a source: " + str(data_vol_least_connected) + " at source: " + str(least_connected_source) + " \nAverage delay: " + str(delay_average_mins) + " mins.  \nDelay by satelite: " + str(delay_allsats_mins))
    file1.close()


end_time = time.time()
runtime = (end_time - start_time)
print("Runtime: %s seconds" % runtime)