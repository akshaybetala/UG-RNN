"""Trains and Evaluates the model."""
# pylint: disable=missing-docstring
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time
import math
from ugrnn import network, utils, input_data
import numpy as np

np.set_printoptions(threshold=np.inf)

from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf
from tensorflow.python.ops import math_ops


import argparse
import matplotlib.pyplot as plt

# Basic model parameters as external flags.
FLAGS = None

class UGRNN(object):
  nn1_hidden_size = [7,7,7,7,7,7,7,7,7,7,3,4,5,6,7,8,9,10,11,12]
  nn1_output_size = [3,4,5,6,7,8,9,10,11,12,3,3,3,3,3,3,3,3,3,3]
  nn2_hidden_size = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
  num_of_networks = 1

  def __init__(self, sess):
    self.feature_pl = tf.placeholder(tf.float32, shape=[None,None, utils.num_of_features()])
    self.path_pl = tf.placeholder(tf.int32, shape=[None,None,3])
    self.target_pl = tf.placeholder(tf.float32)
    self.sequence_len_pl = tf.placeholder(tf.int32)
    self.global_step = tf.Variable(0, name='global_step', trainable=False)
    self.global_step_update_op = tf.assign(self.global_step, tf.add(self.global_step, tf.constant(1)))
    
    # self.learning_rate =  self.linear_rate_decay(FLAGS.max_learning_rate, FLAGS.min_learning_rate) 

    self.learning_rate = tf.train.exponential_decay(learning_rate = FLAGS.learning_rate,
                                                    global_step =  self.global_step, 
                                                    decay_steps = 100, 
                                                    decay_rate = 0.95,
                                                    staircase=True)

    tf.summary.scalar('learning_rate', self.learning_rate)

    self.prediction_ops = []
    self.loss_ops = []
    self.train_ops = []
    self.sess = sess

    for i in xrange(self.num_of_networks):
      network_name = "model_{}".format(i)
      with tf.variable_scope(network_name) as scope:
        # Build a Graph that computes predictions from the inference model.
        model = network.Network(name = network_name,
                                encoding_nn_hidden_size = self.nn1_hidden_size[0],
                                encoding_nn_output_size = self.nn1_output_size[0],
                                output_nn_hidden_size = self.nn2_hidden_size[0],
                                feature_pl = self.feature_pl,
                                path_pl = self.path_pl,
                                target_pl = self.target_pl,
                                sequence_len = self.sequence_len_pl,
                                learning_rate = self.learning_rate,
                                loss_type=FLAGS.loss_type,
                                activation_type=FLAGS.activation_fun,
                                max_seq_len=FLAGS.max_seq_len)
        
        self.prediction_ops.append(model.prediction_op)
        self.loss_ops.append(model.loss_op)
        self.train_ops.append(model.train_op)
  
  def linear_rate_decay(self, max_learning_rate, min_learning_rate):
    step = (max_learning_rate - min_learning_rate) / FLAGS.max_epochs
    return math_ops.sub(max_learning_rate, math_ops.mul(step, tf.to_float(self.global_step)))
    
  '''
  Use the validation set to optimize select the top 10 networks with the minimum RMSE error
  '''
  def optimize(self, dataset):
    n = UGRNN.num_of_networks
    total_loss_values = np.zeros(n)
    while dataset.epochs_completed < 1: 
      feed_dict = self.fill_feed_dict(dataset)
      loss_values, prediction_values = self.sess.run([self.loss_ops, self.prediction_ops],
                                 feed_dict=feed_dict)
      total_loss_values += np.array(loss_values)

    #Get the 10 netowrks with minimum error 
    self.index_of_best_networks = total_loss_values.argsort()[int(-n/2):]
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
                    self.target_pl : targets_feed,
                    self.sequence_len_pl : molecules_feed.feature_vector.shape[1]}
      
      return feed_dict

def main(_):

  
  with tf.Graph().as_default():
    # Create a session for running Ops on the Graph.
    sess = tf.Session()

    print('Creating Graph')
    
    ugrnn_model = UGRNN(sess=sess)
    
    # Build the summary operation based on the TF collection of Summaries.
    summary_op = tf.merge_all_summaries()

    # Create a saver for writing training checkpoints.
    saver = tf.train.Saver()
    
    # Instantiate a SummaryWriter to output summaries and the Graph.
    summary_writer = tf.train.SummaryWriter(FLAGS.train_dir, sess.graph)

    print('Initializing')
    
    # Run the Op to initialize the variables.
    init = tf.initialize_all_variables()
    
    sess.run(init)

    print('Reading Delaney Solubility DataSet')
    
    data_sets = input_data.read_data_sets()
    
    print('Start Training')
    train_dataset = data_sets.train
    validation_dataset = data_sets.validation

    plt.axis([0, FLAGS.max_epochs, 0, 4])
    plt.ion()

    for epoch in xrange(0,FLAGS.max_epochs):
      train_dataset.reset_epoch ()
      total_train_loss = 0
      
      for i in xrange(train_dataset.num_examples):
        feed_dict = ugrnn_model.fill_feed_dict(train_dataset)
        _, l = sess.run([ugrnn_model.train_ops, ugrnn_model.loss_ops],feed_dict=feed_dict)
        total_train_loss += np.mean(l)

      summary = sess.run(summary_op,feed_dict=feed_dict)
      summary_writer.add_summary(summary,epoch)
      summary_writer.flush()
      
      sess.run([ugrnn_model.global_step_update_op])

      train_loss = total_train_loss/data_sets.train.num_examples
      if FLAGS.loss_type is 'rmse':
        train_loss = math.sqrt(train_loss)
      

      
      plt.scatter(epoch, train_loss)
      plt.pause(0.05)

      if epoch % 50 == 0:
        total_validation_loss = 0
        for i in xrange(validation_dataset.num_examples):
          feed_dict = ugrnn_model.fill_feed_dict(validation_dataset)
          l = sess.run([ugrnn_model.loss_ops],feed_dict=feed_dict)
          total_validation_loss += np.mean(l)

        validation_loss = total_validation_loss/validation_dataset.num_examples
        if FLAGS.loss_type is 'rmse':
          validation_loss = math.sqrt(total_validation_loss)
        

        print("Epoch: {:}, Train Loss: {:}, Validation Loss {:}".format(epoch,train_loss, validation_loss))
      else:
        print("Epoch: {:}, Train Loss: {:}".format(epoch, train_loss))
      
    print('Training Finished')
    


    print('Optimize network')
    ugrnn_model.optimize(data_sets.validation)
    predictions = ugrnn_model.predict(data_sets.test)
    error = loss(data_sets.test.labels, predictions)
    print("Loss: {:}".format(error))

if __name__ == '__main__':

  parser = argparse.ArgumentParser()

  parser.add_argument('--fake_data', nargs='?', const=True, type=bool, default=False,
                      help='If true, uses fake data for unit testing.')
  
  parser.add_argument('--max_epochs', type=int, default=4000,
                      help='Number of epochs to run trainer.')
  
  parser.add_argument('--learning_rate', type=float, default=0.00002,
                      help='Initial learning rate')
  
  parser.add_argument('--max_seq_len', type=int, default=100,
                      help='Maximum numbers of atoms, that can be present in a Molecule')
  
  parser.add_argument('--train_dir', type=str, default='train',
                      help='Directory for storing data')
  
  parser.add_argument('--loss_type', type=str, default='rmse',
                      help='The loss function used for training')

  parser.add_argument('--activation_fun', type=str, default='relu',
                      help='Summaries directory')

  parser.add_argument('--dateset', type=str, default='delaney',
                      help='The data used for testing and training')

  parser.add_argument('--summaries_dir', type=str, default='/tmp/ugrnn_logs',
                      help='Summaries directory')

  parser.add_argument('--model', type=int, default=1,
                      help='The model used from the ensemble')

  FLAGS = parser.parse_args()

  tf.app.run()
  