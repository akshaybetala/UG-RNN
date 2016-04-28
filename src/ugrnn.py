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

def inference(rnn_hidden_size,
			rnn_output_size,
			nn_hidden_size,
			feature_placeholder,
			path_placeholder,
      contextual_vector_placeholder,
      sequence_length_placeholder):
  """Build the ugrnn model up to where it may be used for inference."""
  # dims = tf.shape(feature_placeholder).eval()
  # print(dims)
  batch_size = tf.shape(feature_placeholder)[0]

  rnn_input_size = tf.shape(feature_placeholder)[2] + rnn_output_size
  
  rnn_hidden_cell =  tf.nn.rnn_cell.BasicLSTMCell(rnn_hidden_size, forget_bias=0.0, input_size = rnn_input_size)
  rnn_output_cell =  tf.nn.rnn_cell.BasicLSTMCell(rnn_output_size, forget_bias=0.0, input_size =rnn_hidden_size)
  cell = tf.nn.rnn_cell.MultiRNNCell([rnn_hidden_cell,rnn_output_cell])
  states = cell.zero_state(batch_size, tf.float32)

  # for time, input_ in enumerate(inputs):
  #     if time > 0: vs.get_variable_scope().reuse_variables()
  #     # pylint: disable=cell-var-from-loop
  #     def output_state():
  #       return cell(input_, state)
  #     # pylint: enable=cell-var-from-loop
  #     if sequence_length:
  #       (output, state) = control_flow_ops.cond(
  #           time >= max_sequence_length,
  #           lambda: zero_output_state, output_state)
  #     else:
  #       (output, state) = output_state()


  with tf.variable_scope("RNN"):
      for time,x in enumerate(feature_placeholder):
        if time > 0: tf.get_variable_scope().reuse_variables()
        
        t1 = contextual_vector[:,path_placeholder[:,time,0]] 
        t2 = feature_placeholder[:, time, :]
        cell_inputs = tf.concat(1, [t1, t2])
        (cell_output, output_state) = cell(cell_inputs,states)
        contextual_vector[:,path_placeholder[:,time,1]] = cell_output
      rnn_output = tf.reduce_sum(cell_output, 0)

  # Hidden 1
  with tf.name_scope('hidden1'):
    weights = tf.Variable(  
        tf.truncated_normal([rnn_output_size, nn_hidden_size],
                            stddev=1.0 / math.sqrt(float(rnn_output_size))),
        name='weights')
    biases = tf.Variable(tf.zeros([nn_hidden_size]),
                         name='biases')
    hidden1 = tf.tanh(tf.matmul(rnn_output, weights) + biases)
  
  with tf.name_scope('out_linear'):
    weights = tf.Variable(
        tf.truncated_normal([rnn_output, 1],
                            stddev=1.0 / math.sqrt(float(rnn_output))),
        name='weights')
    biases = tf.Variable(tf.zeros([1]),
                         name='biases')
    prediction = tf.matmul(hidden1, weights) + biases
  return predicted


def loss(prediction, target):
  """Calculates the loss from the logits and the labels.

  Returns:
    loss: Loss tensor of type float.
  """
  return tf.nn.l2_loss(prediction-target, name=None)


def training(loss, learning_rate):
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
  tf.scalar_summary(loss.op.name, loss)
  # Create the gradient descent optimizer with the given learning rate.
  optimizer = tf.train.GradientDescentOptimizer(learning_rate)
  # Create a variable to track the global step.
  global_step = tf.Variable(0, name='global_step', trainable=False)
  # Use the optimizer to apply the gradients that minimize the loss
  # (and also increment the global step counter) as a single training step.
  train_op = optimizer.minimize(loss, global_step=global_step)
  return train_op


def evaluation(prediction, target):
  """Evaluate the quality of the logits at predicting the label."""
  return loss(prediction,target)
