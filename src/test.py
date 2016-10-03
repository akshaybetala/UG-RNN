from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import tensorflow as tf
from tensorflow.python.ops import array_ops
import Molecule

smile='CC1OC(CC(O)C1O)OC2C(O)CC(OC2C)OC8C(O)CC(OC7CCC3(C)C(CCC4C3CCC5(C)C(CCC45O)C6=CC(=O)OC6)C7)OC8C'
molecule = Molecule.Molecule(smile)

max_seq_len = 
encoding_nn_output_size = 


with tf.Graph().as_default():
    sess = tf.Session()
    sequence_len = tf.constant(molecule.feature_vector.shape[1])
    contextual_features = tf.range(0, max_seq_len*4)

    flattened_idx_offset = tf.range(0, sequence_len) * max_seq_len *4

    indices = 0 + flattened_idx_offset
    values = [indices,indices,indices,indices]
    indices = tf.pack(values, axis=1, name='pack')
    indices = indices + tf.constant([0,1,2,3])
    indices = tf.reshape(indices,[-1])
    contextual_vector = tf.gather(contextual_features,indices)
    contextual_vector = tf.reshape(contextual_vector , [-1,4*encoding_nn_output_size])

    init = tf.initialize_all_variables()
    sess.run(init)
    print(sess.run([contextual_vector]))