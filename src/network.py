from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time

import numpy as np
import tensorflow as tf
from tensorflow.python.ops import array_ops

"""Builds the UGRNN network.

Implements the inference/loss/training pattern for model building.

1. inference() - Builds the model as far as is required for running the network
forward to make predictions.
2. loss() - Adds to the inference model the layers required to generate loss.
3. training() - Adds to the loss model the Ops required to generate and
apply gradients.

"""

flags = tf.app.flags
FLAGS = flags.FLAGS

class Network(object):
  def __init__(self, encoding_nn_hidden_size, encoding_nn_output_size,
               output_nn_hidden_size, feature_pl, path_pl, sequence_len):
    """Build the ugrnn model up to where it may be used for inference."""
   
    max_seq_len = FLAGS.max_seq_len
    
    flattened_idx_offset = tf.range(0, sequence_len) * max_seq_len * 4
    encoding_nn_input_size = 4*encoding_nn_output_size + FLAGS.initial_feature_vector_size
    
    print('encoding_nn_hidden_size', encoding_nn_hidden_size)
    print('encoding_nn_input_size', encoding_nn_input_size)
    print('encoding_nn_output_size', encoding_nn_output_size)
 
    with tf.variable_scope("EncodingNN") as scope:
      step = tf.constant(0)
      contextual_features = tf.get_variable("contextual_features",
                                            [max_seq_len*max_seq_len*4, encoding_nn_output_size],
                                            dtype=tf.float32,
                                            initializer=tf.constant_initializer(0))
    
      contextual_features = contextual_features.assign(
                                        tf.zeros([max_seq_len*max_seq_len*4, encoding_nn_output_size],
                                          dtype=tf.float32))

      Network.create_variable_for_NN(encoding_nn_input_size, encoding_nn_output_size, encoding_nn_hidden_size)

      _,step,_,_,_,contextual_features,_ = tf.while_loop(Network.cond, Network.body,
                                                      [sequence_len, step, feature_pl, path_pl, 
                                                       flattened_idx_offset, contextual_features, encoding_nn_output_size], 
                                                      parallel_iterations=10, back_prop=True, 
                                                      swap_memory=False, name=None)
      
       # use flattened indices1
      step_contextual_features = Network.get_contextual_feature(contextual_features=contextual_features,
                                                                index=0,
                                                                flattened_idx_offset=flattened_idx_offset,
                                                                encoding_nn_output_size=encoding_nn_output_size )


      begin = tf.get_variable("begin1",[3],dtype=tf.int32)
      begin = tf.scatter_update(begin,1,step,use_locking=None)
      step_feature = tf.squeeze(tf.slice(feature_pl,begin,[-1,1,-1]),squeeze_dims=[1])

      inputs = tf.concat(1,[step_contextual_features,step_feature])
      encodings = Network.single_layer_neural_network(inputs,
                                              encoding_nn_input_size,
                                              encoding_nn_output_size,
                                              encoding_nn_hidden_size)

      molecule_encoding = tf.expand_dims(tf.reduce_sum(tf.tanh(encodings), 0),0)
    
    with tf.variable_scope("OutputNN") as scope:
      self.prediction_op = Network.single_layer_neural_network(molecule_encoding,
                                              encoding_nn_output_size,
                                              1,
                                              output_nn_hidden_size)


  @staticmethod
  def create_variable_for_NN(input_size,output_size,hidden_layer_size):
    with tf.variable_scope('input') as scope:
      weights = tf.get_variable("weights",[input_size, hidden_layer_size],
          initializer=tf.random_normal_initializer())
      
      # Create variable named "biases".
      biases = tf.get_variable("biases", [hidden_layer_size],
          initializer=tf.constant_initializer(0.0))
    
    with tf.variable_scope('hidden') as scope:
      weights = tf.get_variable("weights",[hidden_layer_size,output_size],
          initializer=tf.random_normal_initializer())
      
      # Create variable named "biases".
      biases = tf.get_variable("biases", [output_size],
          initializer=tf.constant_initializer(0.0))
      
  @staticmethod
  def single_layer_neural_network(inputs, 
                                  input_size, 
                                  output_size, 
                                  hidden_layer_size):
      # Create variable named "weights".

    with tf.variable_scope('input') as scope:
      weights = tf.get_variable("weights",[input_size, hidden_layer_size])
      
      # Create variable named "biases".
      biases = tf.get_variable("biases", [hidden_layer_size])

      hidden1 = tf.tanh(tf.matmul(inputs, weights) + biases)
    
    with tf.variable_scope('hidden') as scope:
      weights = tf.get_variable("weights",[hidden_layer_size,output_size])
      
      # Create variable named "biases".
      biases = tf.get_variable("biases", [output_size])
      
      return tf.matmul(hidden1, weights) + biases

  @staticmethod
  def single_layer_neural_network1(inputs):
      # Create variable named "weights".

    with tf.variable_scope('input') as scope:
      weights = tf.get_variable("weights")
      print(tf.shape(inputs))
      # Create variable named "biases".
      biases = tf.get_variable("biases")

      hidden1 = tf.tanh(tf.matmul(inputs, weights) + biases)
    
    with tf.variable_scope('hidden') as scope:
      weights = tf.get_variable("weights")
      
      # Create variable named "biases".
      biases = tf.get_variable("biases")
      
      return tf.matmul(hidden1, weights) + biases

  @staticmethod
  def cond(sequence_len, 
          step,
          feature_pl,
          path_pl,
          flattened_idx_offset,
          contextual_features,
          encoding_nn_output_size):
    return tf.less(step,sequence_len- 1)

  @staticmethod
  def body(sequence_len, 
          step, 
          feature_pl,
          path_pl,
          flattened_idx_offset,
          contextual_features,
          encoding_nn_output_size):
    
    input_begin = tf.get_variable("begin1",[3],dtype=tf.int32,initializer=tf.constant_initializer(0))
    input_begin = tf.scatter_update(input_begin,1,step,use_locking=None)

    step_feature = tf.squeeze(tf.slice(feature_pl,input_begin,[-1,1,-1]))

    input_idx = tf.slice(path_pl, input_begin, [-1,1,1])
    input_idx = tf.reshape(input_idx,[-1])
    max_seq_len = FLAGS.max_seq_len

    output_begin = tf.get_variable("begin2",[3],dtype=tf.int32,initializer=tf.constant_initializer(0))
    output_begin = tf.scatter_update(output_begin,1,step,use_locking=None)
    output_begin = tf.scatter_update(output_begin,2,1,use_locking=None)

    tf.get_variable_scope().reuse_variables()
      
    contextual_features = tf.get_variable("contextual_features")

    step_contextual_features = Network.get_contextual_feature(contextual_features=contextual_features,
                                                              index=input_idx,
                                                              flattened_idx_offset=flattened_idx_offset,
                                                              encoding_nn_output_size=encoding_nn_output_size )
    
    nn_inputs = tf.concat(1,[step_contextual_features,step_feature])
    updated_contextual_vectors = Network.single_layer_neural_network1(nn_inputs)
    updated_contextual_vectors = tf.tanh(updated_contextual_vectors)
    output_idx = tf.squeeze(tf.slice(path_pl, output_begin, [-1,1, 2]))
    
    contextual_features = Network.update_contextual_features(contextual_features=contextual_features,
                                                             indices = output_idx,
                                                             updates = updated_contextual_vectors,
                                                             flattened_idx_offset = flattened_idx_offset)


    with tf.control_dependencies([contextual_features]):
      return (sequence_len, 
              step+1, 
              feature_pl,
              path_pl,
              flattened_idx_offset,
              contextual_features,
              encoding_nn_output_size)

  def add_training_ops(self, learning_rate):
    """Sets up the training Ops.

    Creates a summarizer to track the loss over time in TensorBoard.

    Creates an optimizer and applies the gradients to all trainable variables.

    The Op returned by this function is what must be passed to the
    `sess.run()` call to cause the model to train.

    Args:
      loss: Loss tensor, from loss().
      learning_rate: The learning rate to use for gradient descent.

    Returns:
      train_op: The Op for training.
    """
    # Add a scalar summary for the snapshot loss.

    tf.scalar_summary(self.loss_op.op.name, self.loss_op)
    # Create the gradient descent optimizer with the given learning rate.
    optimizer = tf.train.GradientDescentOptimizer(learning_rate)
    # Create a variable to track the global step.
    global_step = tf.Variable(0, name='global_step', trainable=False)
    # Use the optimizer to apply the gradients that minimize the loss
    # (and also increment the global step counter) as a single training step.
    self.train_op = optimizer.minimize(self.loss_op, global_step=global_step)    

  def add_loss_ops(self, loss_fun, target_pl):

    """Evaluate the quality of the logits at predicting the label."""
    self.loss_op = loss_fun(self.prediction_op, target_pl)
    
  """
  Contextual vector is flatted array
  index is 1D index with
  """
  @staticmethod
  def get_contextual_feature(contextual_features,index,flattened_idx_offset,encoding_nn_output_size):
    indices = index + flattened_idx_offset
    values = [indices,indices,indices,indices]
    indices = tf.pack(values, axis=1, name='pack')
    indices = indices + tf.constant([0,1,2,3])
    indices = tf.reshape(indices,[-1])
    contextual_vector = tf.gather(contextual_features,indices)
    contextual_vector = tf.reshape(contextual_vector , [-1,4*encoding_nn_output_size])
    return contextual_vector

  @staticmethod
  def update_contextual_features(contextual_features,indices,updates,flattened_idx_offset):
    first_indices,second_indices =  tf.split(1, 2, indices)
    indices = tf.squeeze(first_indices + second_indices)
    indices = indices + flattened_idx_offset
    contextual_features = tf.scatter_add(contextual_features,indices,updates,use_locking=None)
    return contextual_features


    