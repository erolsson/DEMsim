import os

import numpy as np


def main():
    directory = os.path.expanduser(r'~/DEMsim/results/porous_electrode/compaction/')
    time = 10
    contact_data = np.genfromtxt(directory + '/contacts/contacts_' + str(time) + '.dou', delimiter=',')
    contact_idx = np.logical_or(contact_data[:, 6] != 0, contact_data[:, 5] > -1e-3)
    contact_data = contact_data[contact_idx, :]
    print(contact_data.shape)


if __name__ == '__main__':
    main()
