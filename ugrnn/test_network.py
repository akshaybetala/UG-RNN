from Molecule import Molecule
import network
import utils

import tensorflow as tf

flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_float('learning_rate', 0.0001, 'Initial learning rate.')
flags.DEFINE_integer('max_epochs',2000,'Number of epochs to run trainer')
flags.DEFINE_string('train_dir', 'train', 'Directory to put the training data.')
flags.DEFINE_integer('initial_feature_vector_size',utils.num_of_features(),'Size of the individual feature for all the nodes' )
flags.DEFINE_integer('max_seq_len',100,'Size of the maximum molecule')


nn1_hidden_size = [7,7,7,7,7,7,7,7,7,7,3,4,5,6,7,8,9,10,11,12]
nn1_output_size = [3,4,5,6,7,8,9,10,11,12,3,3,3,3,3,3,3,3,3,3]
nn2_hidden_size = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
num_of_networks = 1


with tf.Graph().as_default():
	# Create a session for running Ops on the Graph.
	sess = tf.Session()
	feature_pl = tf.placeholder(tf.float32, shape=[None,None,FLAGS.initial_feature_vector_size])
	path_pl = tf.placeholder(tf.int32, shape=[None,None,3])
	targets_pl = tf.placeholder(tf.float32)
	sequence_len_pl = tf.placeholder(tf.int32)

	model = network.Network(name = "model_1",
							encoding_nn_hidden_size = nn1_hidden_size[0],
							encoding_nn_output_size = nn1_output_size[0],
							output_nn_hidden_size = nn2_hidden_size[0],
							feature_pl = feature_pl,
							path_pl = path_pl,
							sequence_len = sequence_len_pl)

	encodings = model.encodings
	inputs = model.inputs
	# Run the Op to initialize the variables.
	init = tf.initialize_all_variables()
	sess.run(init)

	smile = 'C(Cl)(Cl)'
	molecule = Molecule(smile)
	print molecule.feature_vector
	print molecule.directed_graphs
	feed_dict  = {feature_pl : molecule.feature_vector,
				  path_pl : molecule.directed_graphs,
				  targets_pl : 1,
				  sequence_len_pl : molecule.feature_vector.shape[1]}

	enc, inp = sess.run([encodings, inputs],feed_dict=feed_dict)
	print enc
	print inp

