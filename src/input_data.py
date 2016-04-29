from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import gzip
import os
import tempfile

import numpy
from six.moves import urllib
from six.moves import xrange  
import tensorflow as tf
import Molecule
import parse_solubility_data
import math

class DataSet(object):

  def __init__(self, molecules, labels,):

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

  def next_batch(self, batch_size):
    # Return the next `batch_size` examples from this data set.
    start = self._index_in_epoch
    self._index_in_epoch += batch_size
    if self._index_in_epoch > self._num_examples:
      # Finished epoch
      self._epochs_completed += 1
      # Shuffle the data
      perm = numpy.arange(self._num_examples)
      numpy.random.shuffle(perm)
      self._molecules = self._molecules[perm]
      self._labels = self._labels[perm]
      # Start next epoch
      start = 0
      self._index_in_epoch = batch_size
      assert batch_size <= self._num_examples
    end = self._index_in_epoch
    return self._molecules[start:end], self._labels[start:end]

def extract_molecules_from_smiles(SMILES):
    molecules = []
    size = len(SMILES)
    # TODO: Change it back to batch_size
    for i in xrange(2):
      m = Molecule.Molecule(SMILES[i])
      molecules.append(m)
    return molecules

def read_data_sets():
  class DataSets(object):
    pass
  data_sets = DataSets()

  SMILES, prediction_targets = parse_solubility_data.load_solubility_data(path = '../data/Delaney_solubility.txt')
  molecules = extract_molecules_from_smiles(SMILES)
  size = len(SMILES)
  
  TRAIN_SIZE = int(math.floor(0.7*size))
  VALIDATION_SIZE = int(math.floor(0.1*size))
  train_molecules = molecules[:TRAIN_SIZE]
  train_labels = prediction_targets[:TRAIN_SIZE]

  validation_molecules = molecules[TRAIN_SIZE:TRAIN_SIZE+VALIDATION_SIZE]
  validation_labels = prediction_targets[TRAIN_SIZE:TRAIN_SIZE+VALIDATION_SIZE]

  test_molecules = molecules[TRAIN_SIZE+VALIDATION_SIZE:]
  test_labels = prediction_targets[TRAIN_SIZE+VALIDATION_SIZE:]

  data_sets.train = DataSet(train_molecules, train_labels)
  data_sets.validation = DataSet(validation_molecules, validation_labels)
  data_sets.test = DataSet(test_molecules, test_labels)

  return data_sets

