from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import gzip
import os
import tempfile

import numpy
from six.moves import urllib
from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf
import Molecule

class DataSet(object):

  def __init__(self, molecules, labels,):

    self._num_examples = len(molecules)

    self._molecules = molecules
    self._labels = labels
    self._epochs_completed = 0
    self._index_in_epoch = 0

  @property
  def images(self):
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

  def next_batch(self, batch_size, fake_data=False):
    """Return the next `batch_size` examples from this data set."""
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
    for i in xrange(size):
      molecules.append(Molecule(SMILES[i]))
    return molecules

def read_data_sets(dtype=tf.float32):
  class DataSets(object):
    pass
  data_sets = DataSets()

  SMILES, prediction_targets = load_solubility_data(path = '../data/Delaney_solubility.txt')
  molecules = extract_molecules_from_smiles(SMILES)
  size = len(SMILES)
  
  TRAIN_SIZE = 0.7*size
  VALIDATION_SIZE = 0.1*size
  train_molecules = molecules[:TRAIN_SIZE]
  train_labels = prediction_targets[:TRAIN_SIZE]

  validation_molecules = molecules[TRAIN_SIZE:TRAIN_SIZE+VALIDATION_SIZE]
  validation_labels = prediction_targets[TRAIN_SIZE:TRAIN_SIZE+VALIDATION_SIZE]

  test_molecules = molecules[TRAIN_SIZE+VALIDATION_SIZE:]
  test_labels = prediction_targets[TRAIN_SIZE+VALIDATION_SIZE:]

  data_sets.train = DataSet(train_images, train_labels, dtype=dtype)
  data_sets.validation = DataSet(validation_images, validation_labels,
                                 dtype=dtype)
  data_sets.test = DataSet(test_images, test_labels, dtype=dtype)

  return data_sets