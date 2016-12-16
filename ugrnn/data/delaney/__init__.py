from __future__ import print_function

import os

import numpy as np


def read_data_set():
    path = os.path.dirname(os.path.realpath(__file__))
    smiles_path = os.path.join(path, 'delaney.smi')
    target_path = os.path.join(path, 'delaney.target')
    perm_smiles_path = os.path.join(path, 'delaney_perm.smi')
    perm_target_path = os.path.join(path, 'delaney_perm.target')

    with open(perm_smiles_path) as f:
        smiles = np.array([line.rstrip() for line in f])

    with open(perm_target_path) as f:
        labels = np.array([float(line.rstrip()) for line in f])

    return smiles, labels
