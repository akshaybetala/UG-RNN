from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from six.moves import xrange

from ugrnn.molecule import Molecule
from ugrnn.utils import read_csv, permute_data
import csv

class DataSet(object):

    def __init__(self, csv_file_path, smile_col_name, target_col_name, logp_col_name=None, contract_rings=False):
        data = read_csv(csv_file_path, smile_col_name, target_col_name, logp_col_name)

        smiles = np.array(zip(*data)[0])
        self._labels = np.array(zip(*data)[1])
        logp = np.array(zip(*data)[2]) if logp_col_name else None

        self._num_examples = len(smiles)
        self._molecules = extract_molecules_from_smiles(smiles, logp, contract_rings)
        self._epochs_completed = 0
        self._index_in_epoch = 0

    @property
    def molecules(self):
        return self._molecules

    @property
    def labels(self):
        return self._labels

    @property
    def num_examples(self):
        return self._num_examples

    @property
    def epochs_completed(self):
        return self._epochs_completed

    @property
    def index_in_epoch(self):
        return self._index_in_epoch

    def reset_epoch(self, permute=False):
        self._index_in_epoch = 0
        self._epochs_completed = 0

        if permute:
            self.permute_data()


    def next_batch(self, batch_size):
        """Return the next `batch_size` examples from this data set."""

        start = self._index_in_epoch
        self._index_in_epoch += batch_size
        if self._index_in_epoch > self._num_examples:
            # Finished epoch
            self._epochs_completed += 1
            # Shuffle the data
            self.permute_data()
            # Start next epoch
            start = 0
            self._index_in_epoch = batch_size
            assert batch_size <= self._num_examples
        end = self._index_in_epoch
        return self._molecules[start:end], self._labels[start:end]

    def permute_data(self):
        perm = np.random.permutation(self._num_examples)
        self._molecules = self._molecules[perm]
        self._labels = self._labels[perm]


def extract_molecules_from_smiles(SMILES, logp, contract_rings):
    size = len(SMILES)
    molecules = np.empty(size, dtype=object)
    for i in xrange(size):
        molecule_logp = None
        if logp is not None:
            molecule_logp = logp[i]
        molecules[i] = Molecule(SMILES[i], molecule_logp, contract_rings)
    return molecules
