"""Trains and Evaluates the model."""
# pylint: disable=missing-docstring
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time
from ugrnn import network, utils, input_data
import numpy as np

np.set_printoptions(threshold=np.inf)

from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf

# Basic model parameters as external flags.
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_float('learning_rate', 0.0001, 'Initial learning rate.')
flags.DEFINE_integer('max_epochs',2000,'Number of epochs to run trainer')
flags.DEFINE_string('train_dir', 'train', 'Directory to put the training data.')
flags.DEFINE_integer('initial_feature_vector_size',utils.num_of_features(),'Size of the individual feature for all the nodes' )
flags.DEFINE_integer('max_seq_len',100,'Size of the maximum molecule')

class UGRNN(object):
  nn1_hidden_size = [7,7,7,7,7,7,7,7,7,7,3,4,5,6,7,8,9,10,11,12]
  nn1_output_size = [3,4,5,6,7,8,9,10,11,12,3,3,3,3,3,3,3,3,3,3]
  nn2_hidden_size = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
  num_of_networks = 20

  def __init__(self, sess, global_step, loss_fun):
    self.feature_pl = tf.placeholder(tf.float32, shape=[None,None,FLAGS.initial_feature_vector_size])
    self.path_pl = tf.placeholder(tf.int32, shape=[None,None,3])
    self.targets_pl = tf.placeholder(tf.float32)
    self.sequence_len_pl = tf.placeholder(tf.int32)
    self.global_step = global_step

    self.prediction_ops = []
    self.loss_ops = []
    self.train_ops = []
    self.sess = sess

    for i in xrange(self.num_of_networks):
      with tf.variable_scope("model_{}".format(i)) as scope:
        # Build a Graph that computes predictions from the inference model.
        model = network.Network(encoding_nn_hidden_size = self.nn1_hidden_size[i],
                                encoding_nn_output_size = self.nn1_output_size[i],
                                output_nn_hidden_size = self.nn2_hidden_size[i],
                                feature_pl = self.feature_pl,
                                path_pl = self.path_pl,
                                sequence_len = self.sequence_len_pl)
        
        model.add_loss_ops(loss_fun=loss_fun, target_pl = self.targets_pl)
        model.add_training_ops(global_step=global_step, initial_learning_rate=FLAGS.learning_rate)
        self.prediction_ops.append(model.prediction_op)
        self.loss_ops.append(model.loss_op)
        self.train_ops.append(model.train_op)
    
  def train(self, dataset, epochs=1):
    # start the training loop.

    for i in xrange(epochs):
      dataset.reset_epoch ()
      while dataset.epochs_completed < 1:
        feed_dict = self.fill_feed_dict(dataset)
        _ = self.sess.run([self.train_ops],feed_dict=feed_dict)
      update = tf.assign(self.global_step, tf.add(self.global_step, tf.constant(1)))
      self.sess.run([update])

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

def aae(predictions, targets):
  return np.abs(predictions - targets).mean()

def rmse_loss(predictions, targets):
  # return tf.nn.l2_loss(prediction-targets, name=None)
  return tf.reduce_mean(tf.square(predictions-targets))

def aae_loss(predictions, targets):
  return tf.reduce_mean(tf.abs(predictions-targets))

def main(_):

  print('Reading Delaney Solubility DataSet')
  data_sets = input_data.read_data_sets()
    

  loss_type = 'rmse'
  
  if loss_type is 'aae':
    loss_fun = aae_loss
    loss = aae
  else:
    loss_fun = rmse_loss
    loss = rmse
  
  with tf.Graph().as_default():
    # Create a session for running Ops on the Graph.
    sess = tf.Session()

    global_step = tf.Variable(0, name='global_step', trainable=False)

    print('Creating Graph')
    ugrnn_model = UGRNN(sess=sess, global_step=global_step, loss_fun=loss_fun)
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

    print('Start Training')
    EPOCHS =0
    epochs_per_train = 5
    while EPOCHS < FLAGS.max_epochs:
      ugrnn_model.train(dataset=data_sets.train,epochs=epochs_per_train)
      EPOCHS+=epochs_per_train
      
      predictions = ugrnn_model.predict(data_sets.train)
      fp = open("results",'w')
      fp.write('\n'.join('%s %s' % x for x in zip(data_sets.train.labels, predictions)))
      fp.close()
      error1 = loss(data_sets.train.labels, predictions)

      predictions = ugrnn_model.predict(data_sets.validation)
      error2 = loss(data_sets.validation.labels, predictions)
      
      print("Epoch: {:}, Train Loss: {:}, Validation Loss: {:}".format(EPOCHS,error1,error2))

    print('Training Finished')
    print('Optimize network')
    ugrnn_model.optimize(data_sets.validation)
    predictions = ugrnn_model.predict(data_sets.test)
    error = loss(data_sets.test.labels, predictions)
    print("Loss: {:}".format(error))
if __name__ == '__main__':
  tf.app.run()
  