import os

import numpy as np


def read_data_set():
    path = os.path.dirname(os.path.realpath(__file__))
    smiles_path = os.path.join(path, 'delaney.smi')
    target_path = os.path.join(path, 'delaney.target')

    with open(smiles_path) as f:
        smiles = np.array([line.rstrip() for line in f])

    with open(target_path) as f:
        labels = np.array([float(line.rstrip()) for line in f])

    return smiles, labels
