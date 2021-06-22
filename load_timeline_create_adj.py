import numpy as np
import pandas as pd
import time

start_time = time.time()
case_name = "random_test"
timeline_file_location = "Timeline_files\\"
adj_file_location = "Adjacency_matrices\\"


# define number of sources, satellites and sinks
num_sources = 300                                           
num_sats = 100
num_gstations = 77

# number of days simulated
num_days = 100

#set timestep length in seconds
time_step_length = 30

#number of timesteps to consider
num_timesteps = int(num_days * 24 * 60 * 60 / time_step_length)

# matrix size is number of sources + satellites + sinks along rows and columns
num_matrix_elements = num_sources + num_sats + num_gstations

adjacency_matrix_basic = np.zeros((num_matrix_elements, num_matrix_elements))

# to create adjacency matrix
for sat_number in range(0,num_sats):
    print("adj_matrix_" + str(sat_number))
    timeline = np.load(timeline_file_location + 'timeline_sat_' + str(sat_number) + '_' + case_name + '.npy', allow_pickle = True)
    timestep_number = 0                                      # set timestep to count from 0 for new satellite
    for step in timeline:                                    # for each step in the timeline
        for element in step[1:]:                                    # for each item seen in this time step
            if element[0:8] == 'Gstation':                                                   
                gstation_number_temp = int(element[9:])          # pull out the ground station number
                # add the time step length to the weighting at the element row = satellite, column = ground station
                adjacency_matrix_basic[num_sources + sat_number, num_sources + num_sats + gstation_number_temp - 1] += time_step_length
            elif element[0:6] == 'Source':                          # if it's a target 
                target_number_temp = int(element[7:])          # pull out the target number
                # add the time step length to the weighting at the element row = target, column = satellite
                adjacency_matrix_basic[target_number_temp -1, num_sources + sat_number] += time_step_length


        timestep_number += 1                                        # increment timestep

adjacency_matrix_basic_norm = adjacency_matrix_basic / num_days

# save adjacency matrices
data_frame_temp = pd.DataFrame(adjacency_matrix_basic)                          # convert to data frame using pandas
data_frame_temp.to_csv(adj_file_location + 'adj_matrix_basic_' + case_name + '_' + str(num_days) + '.csv', index=False, header= None)   # save as csv file          
data_frame_temp = pd.DataFrame(adjacency_matrix_basic_norm)                          # convert to data frame using pandas
data_frame_temp.to_csv(adj_file_location + 'adj_matrix_basic_norm' + case_name + '_' + str(num_days) + '.csv', index=False, header= None)   # save as csv file          


end_time = time.time()
runtime = (end_time - start_time)
print("Runtime: %s seconds" % runtime)