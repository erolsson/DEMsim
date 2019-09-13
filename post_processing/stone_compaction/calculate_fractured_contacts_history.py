import pickle
import sys

import numpy as np

base_directory = sys.argv[-1]
surface_position_file = base_directory + "/compaction/surface_positions.dou"
surface_forces_file = base_directory + "/compaction/surface_forces.dou"

position_data = np.genfromtxt(surface_position_file, delimiter=",", usecols=(-1, -2))
data = np.zeros((position_data.shape[0], 4))
data[:, 0:2] = position_data

force_data = np.genfromtxt(surface_forces_file, delimiter=",", usecols=(-2,))
data[:, 2] = force_data
start_idx = -1
current_time = 0.
active_cracks = set()
idx = -1
for line in open(base_directory + "/compaction/particle_cracks.dou"):
    crack_data = line.split(",")
    time = float(crack_data[-1])
    if idx == -1:
        idx = np.where(data[:, 0] == time)[0][0] - 1  # -1 as the index is increased in the next if statement

    if time != current_time:
        if time.is_integer():
            print time
        active_cracks = set()
        current_time = time
        idx += 1

    crack_info = (int(crack_data[0][3:]), int(crack_data[1][13:]), float(crack_data[-4]),
                  float(crack_data[-3]), float(crack_data[-2]))
    if crack_info not in active_cracks:
        data[idx, 3] += 1
        active_cracks.add(crack_info)

with open(base_directory + "/crack_history.pkl", 'w') as pickle_handle:
    pickle.dump(data, pickle_handle)
