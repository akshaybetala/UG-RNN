from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from six.moves import xrange

from ugrnn.data import delaney
from ugrnn.molecule import Molecule


class DataSet(object):
    def __init__(self, molecules, labels, ):
        self._num_examples = len(molecules)

        self._molecules = molecules
        self._labels = labels
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
            self._molecules, self._labels = permute_data(self._molecules, self._labels)

    def next_batch(self, batch_size):
        """Return the next `batch_size` examples from this data set."""

        start = self._index_in_epoch
        self._index_in_epoch += batch_size
        if self._index_in_epoch > self._num_examples:
            # Finished epoch
            self._epochs_completed += 1
            # Shuffle the data
            self._molecules, self._labels = permute_data(self._molecules, self._labels)
            # Start next epoch
            start = 0
            self._index_in_epoch = batch_size
            assert batch_size <= self._num_examples
        end = self._index_in_epoch
        return self._molecules[start:end], self._labels[start:end]


def extract_molecules_from_smiles(SMILES, contract_rings):
    size = len(SMILES)
    molecules = np.empty(size, dtype=object)
    for i in xrange(size):
        molecules[i] = Molecule(SMILES[i], contract_rings)
    return molecules


def permute_data(data, labels):
    data_len = len(data)
    perm = np.random.permutation(data_len)
    data_perm = data[perm]
    labels_perm = labels[perm]
    return data_perm, labels_perm


def cross_validation_split(data, labels, crossval_split_index, crossval_total_num_splits, validation_data_ratio=0.1):
    '''
    Manages cross-validation splits given fixed lists of data/labels
    <crossval_total_num_splits> directly affects the size of the test set ( it is <size of data-set>/crossval_total_num_splits)
    Returns:
    ----------
        traindata, valdata, testdata

    '''

    assert validation_data_ratio > 0 and validation_data_ratio < 1
    assert crossval_split_index < crossval_total_num_splits

    N = len(data)
    n_test = int(N * 1. / crossval_total_num_splits)
    if crossval_split_index == crossval_total_num_splits - 1:
        n_test = N - crossval_split_index * n_test

    start_test = crossval_split_index * n_test
    end_test = crossval_split_index * n_test + n_test
    testdata = (data[start_test: end_test], labels[start_test: end_test])

    rest_data = np.concatenate((data[:start_test], data[end_test:]), axis=0)
    rest_labels = np.concatenate((labels[:start_test], labels[end_test:]), axis=0)

    n_valid = int(N * validation_data_ratio)
    valdata = (rest_data[: n_valid], rest_labels[: n_valid])
    traindata = (rest_data[n_valid:], rest_labels[n_valid:])
    print(len(traindata[0]), len(valdata[0]), len(testdata[0]))
    return traindata, valdata, testdata

def read_data_sets(dataset="delaney", contract_rings=False):
    class DataSets(object):
        pass

    data_sets = DataSets()
    if dataset == "delaney":
        smiles, labels = delaney.read_data_set()

    molecules = extract_molecules_from_smiles(smiles, contract_rings)

    traindata, valdata, testdata = cross_validation_split(data=molecules,
                                                          labels=labels,
                                                          crossval_split_index=0,
                                                          crossval_total_num_splits=10,
                                                          validation_data_ratio=0.1)

    data_sets.train = DataSet(traindata[0], traindata[1])
    data_sets.validation = DataSet(valdata[0], valdata[1])
    data_sets.test = DataSet(testdata[0], testdata[1])

    return data_sets
