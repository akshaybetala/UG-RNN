"""Trains and Evaluates the model."""
# pylint: disable=missing-docstring
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time
import input_data
import network
import utils
import numpy as np
np.set_printoptions(threshold=np.inf)

from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf


# Basic model parameters as external flags.
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_float('learning_rate', 0.01, 'Initial learning rate.')
flags.DEFINE_integer('max_epochs',5000,'Number of epochs to run trainer')
flags.DEFINE_string('train_dir', 'train', 'Directory to put the training data.')
flags.DEFINE_integer('initial_feature_vector_size',utils.num_of_features(),'Size of the individual feature for all the nodes' )
flags.DEFINE_integer('max_seq_len',100,'Size of the maximum molecule')


class UGRNN(object):
  nn1_hidden_size = [7,7,7,7,7,7,7,7,7,7,3,4,5,6,7,8,9,10,11,12] 
  nn1_output_size = [3,4,5,6,7,8,9,10,11,12,3,3,3,3,3,3,3,3,3,3]
  nn2_hidden_size = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
  num_of_networks = 20

  def __init__(self, sess):
    self.feature_pl = tf.placeholder(tf.float32, shape=[None,None,FLAGS.initial_feature_vector_size])
    self.path_pl = tf.placeholder(tf.int32, shape=[None,None,3])
    self.targets_pl = tf.placeholder(tf.float32)
    self.sequence_len_pl = tf.placeholder(tf.int32)
    
    self.prediction_ops = []
    self.loss_ops = []
    self.train_ops = []

    self.sess = sess

    for i in xrange(20):
      with tf.variable_scope("model_{}".format(i)) as scope:
        # Build a Graph that computes predictions from the inference model.
        model = network.Network(encoding_nn_hidden_size = self.nn1_hidden_size[i],
                                encoding_nn_output_size = self.nn1_output_size[i],
                                output_nn_hidden_size = self.nn2_hidden_size[i],
                                feature_pl = self.feature_pl,
                                path_pl = self.path_pl,
                                sequence_len = self.sequence_len_pl)
        
        model.add_loss_ops(loss_fun=rmse_loss, target_pl = self.targets_pl)
        model.add_training_ops(FLAGS.learning_rate)

        self.prediction_ops.append(model.prediction_op)
        self.loss_ops.append(model.loss_op)
        self.train_ops.append(model.train_op)



  def train(self, dataset, epochs=1):
    # start the training loop.
    dataset.reset_epoch ()
    while dataset.epochs_completed < epochs:

      # Fill a feed dictionary with the molecule and target
      # for this particular training step.
      feed_dict = self.fill_feed_dict(dataset)
      _, loss_values, prediction_values = self.sess.run([self.train_ops, self.loss_ops, self.prediction_ops],
                               feed_dict=feed_dict)
        
  '''
  Use the validation set to optimize select the top 10 networks with the minimum RMSE error
  '''
  def optimize(self, dataset):
    total_loss_values = np.zeros(20)
    while dataset.epochs_completed < 1:
      feed_dict = self.fill_feed_dict(dataset)
      loss_values, prediction_values = self.sess.run([self.loss_ops, self.prediction_ops],
                                 feed_dict=feed_dict)
      total_loss_values += np.array(loss_values)

    #Get the 10 netowrks with minimum error 
    self.index_of_best_networks = total_loss_values.argsort()[:10]
    self.final_prediction_ops = [ self.prediction_ops[index] for index in self.index_of_best_networks]
  

  def predict(self, dataset):
    dataset.reset_epoch()
    predictions = []
    while dataset.epochs_completed < 1:
      feed_dict = self.fill_feed_dict(dataset)
      prediction_values = self.sess.run([self.prediction_ops], feed_dict=feed_dict)
      predictions.append(np.mean(prediction_values))
    
    return np.array(predictions)
        
  def fill_feed_dict(self, dataset):
      molecules_feed, targets_feed = dataset.next_molecule()
      feed_dict  = {self.feature_pl : molecules_feed.feature_vector,
                self.path_pl : molecules_feed.directed_graphs,
                self.targets_pl : targets_feed,
                self.sequence_len_pl : molecules_feed.feature_vector.shape[1]}
                
      return feed_dict

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def rmse_loss(prediction, target):
    """Calculates the loss from the logits and the labels.

    Returns:
      loss: Loss tensor of type float.
    """
    return tf.nn.l2_loss(prediction-target, name=None)

def main(_):

  
  with tf.Graph().as_default():
    # Create a session for running Ops on the Graph.
    sess = tf.Session()

    ugrnn_model = UGRNN(sess)
    # Build the summary operation based on the TF collection of Summaries.
    summary_op = tf.merge_all_summaries()

    # Create a saver for writing training checkpoints.
    saver = tf.train.Saver()
    
    # Instantiate a SummaryWriter to output summaries and the Graph.
    summary_writer = tf.train.SummaryWriter(FLAGS.train_dir, sess.graph)

    # Run the Op to initialize the variables.
    init = tf.initialize_all_variables()
    sess.run(init)

    data_sets = input_data.read_data_sets()
    EPOCHS =0
    while EPOCHS < FLAGS.max_epochs:
      ugrnn_model.train(dataset=data_sets.train,epochs=2)
      EPOCHS+=2
      predictions = ugrnn_model.predict(data_sets.test)
      error = rmse(data_sets.test.labels, predictions)
      print("Epoch: {:}, Loss: {:}".format(EPOCHS,error))

    ugrnn_model.optimize(data_sets.validation)
    predictions = ugrnn_model.predict(data_sets.test)
    error = rmse(data_sets.test.labels, predictions)

if __name__ == '__main__':
  tf.app.run()
