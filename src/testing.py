from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time

from six.moves import xrange
import tensorflow as tf

import input_data
import ugrnn
import utils


input_data = np.ones((1000,100,rnn_input_size))
sequence_data = np.ones()
batch_size = 100


rnn_input_size = 10

x = tf.placeholder(tf.float32, shape=(None,None,rnn_input_size)) 
sequence_lengths = tf.placeholder(tf.int32,shape = (None))

lstm_cell =  tf.nn.rnn_cell.BasicLSTMCell(rnn_input_size, forget_bias=0.0)
state = lstm_cell.zero_state(batch_size,tf.float32)
outputs, final_state = tf.nn.rnn(lstm_cell, x, state, sequence_length = sequence_lengths)

# code for one epoch
iterations = total_data_length / batch_size
max_sequence_length = max(all_possible_sequence_lengths)
cur_state = initial_state

for i in range(iterations):
    # x is of dimension [max_sequence_length, batch_size, input_size]
    # sequence_lengths is of dimension [batch_size]
    x_data, sequence_data, y_data = mini_batch(batch_size)

    feed_dict = {k: v for k, v in zip(x, x_data)}
    feed_dict.append(sequence_lengths: sequence_data, ...)
    outs, cur_state, _ = session.run([outputs, final_state, train], feed_dict)