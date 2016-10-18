from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import math
from six.moves import xrange  
import tensorflow as tf
from ugrnn.Molecule import Molecule
from ugrnn.data import delaney

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

  @property
  def index_in_epoch(self):
    return self._index_in_epoch

  def reset_epoch(self):
    self._index_in_epoch = 0
    self._epochs_completed = 0

    perm = np.arange(self._num_examples)
    np.random.shuffle(perm)
    self._molecules = self._molecules[perm]
    self._labels = self._labels[perm]

  def next_molecule(self):
    # Return the next example from this data set.
    
    # Finished epoch
    if self._index_in_epoch == self._num_examples-1:
      self._epochs_completed += 1
      # Start next epoch
      self._index_in_epoch = 0
      return self._molecules[self._num_examples-1], self._labels[self._num_examples-1]
    self._index_in_epoch +=1
    return self._molecules[self._index_in_epoch-1], self._labels[self._index_in_epoch-1]

def extract_molecules_from_smiles(SMILES):
    size = len(SMILES)
    molecules = np.empty(size,dtype=object)
    for i in xrange(size):
      molecules[i] = Molecule(SMILES[i])
    return molecules

def read_data_sets(dataset_type = 'delaney'):

  class DataSets(object):
    pass
  data_sets = DataSets()

  smiles_path = delaney.smiles_path
  target_path = delaney.target_path

  with open(smiles_path) as f:
    smiles = np.array([line.rstrip() for line in f])

  with open(target_path) as f:
    prediction_targets = np.array([float(line.rstrip()) for line in f])
    

  molecules = extract_molecules_from_smiles(smiles)
  num_examples = len(smiles)

  total_training_size = int(math.floor(0.9 * num_examples) )

  total_train_molecules = molecules[:total_training_size]
  total_train_labels = prediction_targets[:total_training_size]

  test_molecules = molecules[total_training_size:]
  test_labels = prediction_targets[total_training_size:]

  #Now divide the total training in training and validation set
  if num_examples > 100:
    train_ratio = .8
  else:
    train_ratio = .85

  train_size = int(train_ratio*total_training_size)

  train_molecules = total_train_molecules[:train_size]
  train_labels = total_train_labels[:train_size]

  validation_molecules = total_train_molecules[train_size:]
  validation_labels = total_train_labels[train_size:]

  data_sets.train = DataSet(train_molecules, train_labels)
  data_sets.validation = DataSet(validation_molecules, validation_labels)
  data_sets.test = DataSet(test_molecules, test_labels)

  return data_sets

