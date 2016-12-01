import tensorflow as tf

import network
import logging
from molecule import Molecule

flags = tf.app.flags
flags.DEFINE_integer('initial_feature_vector_size', Molecule.num_of_features(),
                     'Size of the individual feature for all the nodes')
flags.DEFINE_integer('max_seq_len', 100, 'Size of the maximum molecule')
flags.DEFINE_string('loss_type', 'rmse','')
flags.DEFINE_string('activation_type', 'identity','')
FLAGS = flags.FLAGS

log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_format)
logger = logging.getLogger(__name__)

network.FLAGS = FLAGS

with tf.Graph().as_default():
    # Create a session for running Ops on the Graph.
    sess = tf.Session()
    feature_pl = tf.placeholder(tf.float32, shape=[None, None, FLAGS.initial_feature_vector_size])
    path_pl = tf.placeholder(tf.int32, shape=[None, None, 3])
    target_pl = tf.placeholder(tf.float32)
    sequence_len_pl = tf.placeholder(tf.int32)

    model = network.Network(name="test_model",
                            encoding_nn_hidden_size=2,
                            encoding_nn_output_size=1,
                            output_nn_hidden_size=1,
                            feature_pl=feature_pl,
                            path_pl=path_pl,
                            target_pl=target_pl,
                            sequence_len=sequence_len_pl,
                            learning_rate=0.0001,
                            max_seq_len=100,
                            initializer="one")

    prediction = model.prediction_op

    # Run the Op to initialize the variables.
    init = tf.global_variables_initializer()
    sess.run(init)

    smile = 'C(Cl)'
    molecule = Molecule(smile)
    print(molecule.directed_graphs)
    print(molecule.feature_vector)
    feed_dict = {feature_pl: molecule.feature_vector,
                 path_pl: molecule.directed_graphs,
                 target_pl: 1,
                 sequence_len_pl: molecule.feature_vector.shape[1]}

    p = sess.run([prediction], feed_dict=feed_dict)
    print(p)
