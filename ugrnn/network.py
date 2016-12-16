from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import logging
import tensorflow as tf
from ugrnn.molecule import Molecule
import os
logger = logging.getLogger(__name__)

FLAGS = None


class Network(object):
    max_seq_len=0
    def __init__(self, name, encoding_nn_hidden_size, encoding_nn_output_size,
                 output_nn_hidden_size, max_seq_len, feature_pls, path_pls, sequence_lens,
                 target_pls, learning_rate, initializer):
        """Build the ugrnn model up to where it may be used for inference."""
        self.name = name
        self.encoding_nn_output_size = encoding_nn_output_size
        self.output_nn_hidden_size = output_nn_hidden_size
        Network.max_seq_len = max_seq_len
        self.encoding_nn_hidden_size = encoding_nn_hidden_size

        self.encoding_nn_input_size = 4 * encoding_nn_output_size + Molecule.num_of_features()
        self.learning_rate = learning_rate
        self.initializer_fun = Network.get_initializer(initializer)

        assert FLAGS.batch_size == len(feature_pls)
        assert FLAGS.batch_size == len(path_pls)
        assert FLAGS.batch_size == len(target_pls)
        assert FLAGS.batch_size == len(sequence_lens)

        self.feature_pls = feature_pls
        self.path_pls = path_pls
        self.target_pls = target_pls
        self.sequence_lens = sequence_lens
        self.flattened_idx_offsets = [(tf.range(0, sequence_lens[i]) * Network.max_seq_len * 4) for i in
                                      xrange(0, FLAGS.batch_size)]

        self.trainable_variables = []
        self.create_network_variable()

        self.prediction_ops = [self.predict(self.feature_pls[0], self.path_pls[0],
                                            self.sequence_lens[0], self.flattened_idx_offsets[0])]
        for i in xrange(1, FLAGS.batch_size):
            with tf.control_dependencies([self.prediction_ops[i - 1]]):
                prediction_op = self.predict(self.feature_pls[i], self.path_pls[i],
                                             self.sequence_lens[i], self.flattened_idx_offsets[i])
                self.prediction_ops.append(prediction_op)

        self.loss()
        self.training()

    def save_network(self, sess, step):
        self.saver = tf.train.Saver(self.trainable_variables, max_to_keep=1)
        checkpoint_dir = self.get_checkpoint_dir()
        file_path = os.path.join(checkpoint_dir, "checkpoint")
        self.saver.save(sess, save_path=file_path,global_step=step,)

    def restore_network(self, sess):
        saver = tf.train.Saver(self.trainable_variables)
        checkpoint_dir = self.get_checkpoint_dir()
        saver.restore(sess, tf.train.latest_checkpoint(checkpoint_dir))

    def get_checkpoint_dir(self):
        '''
        :return:
        '''
        path = os.path.dirname(os.path.realpath(__file__))
        input_type = "ugrnn" if not FLAGS.contract_rings else "ugrnn_cr"
        folder = "{:}/{:}/{:}/".format(FLAGS.models_dir, self.name, input_type)
        folder_path = os.path.join(path, folder)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        return folder_path

    def create_network_variable(self):
        '''
        :return:
        '''
        with tf.variable_scope("EncodingNN") as scope:
            contextual_features = tf.get_variable("contextual_features",
                                                  [self.max_seq_len * self.max_seq_len * 4,
                                                   self.encoding_nn_output_size],
                                                  dtype=tf.float32,
                                                  initializer=tf.constant_initializer(0),
                                                  trainable=False)

            with tf.variable_scope('hidden1') as scope:
                weights = Network.weight_variable([self.encoding_nn_input_size, self.encoding_nn_hidden_size],
                                                  initializer=self.initializer_fun)
                biases = Network.bias_variable([self.encoding_nn_hidden_size])
                self.trainable_variables.append(weights)
                self.trainable_variables.append(biases)

            with tf.variable_scope('output') as scope:
                weights = Network.weight_variable([self.encoding_nn_hidden_size, self.encoding_nn_output_size],
                                                  initializer=self.initializer_fun)

                biases = Network.bias_variable([self.encoding_nn_output_size])
                self.trainable_variables.append(weights)
                self.trainable_variables.append(biases)

        with tf.variable_scope("OutputNN") as scope:
            with tf.variable_scope('hidden1') as scope:
                weights = Network.weight_variable([self.encoding_nn_output_size, self.output_nn_hidden_size],
                                                  self.initializer_fun, 'weights_decay')

                biases = Network.bias_variable([self.output_nn_hidden_size])
                self.trainable_variables.append(weights)
                self.trainable_variables.append(biases)

            with tf.variable_scope('output') as scope:
                weights = Network.weight_variable([self.output_nn_hidden_size, 1],
                                                  self.initializer_fun, 'weights_decay')
                self.trainable_variables.append(weights)

    def predict(self, feature_pl, path_pl, sequence_len,flattened_idx_offset):
        '''
        :return:
        '''
        with tf.variable_scope("EncodingNN", reuse=True) as scope:
            step = tf.constant(0)
            contextual_features = tf.get_variable("contextual_features")
            contextual_features = contextual_features.assign(
                tf.zeros([self.max_seq_len * self.max_seq_len * 4, self.encoding_nn_output_size],
                         dtype=tf.float32))

            _, step, _, _, _, contextual_features, _ = tf.while_loop(Network.cond, Network.body,
                                                                     [sequence_len, step, feature_pl,
                                                                      path_pl,
                                                                      flattened_idx_offset, contextual_features,
                                                                      self.encoding_nn_output_size],
                                                                     back_prop=True,
                                                                     swap_memory=False, name=None)

            # use flattened indices1
            step_contextual_features = Network.get_contextual_feature(contextual_features=contextual_features,
                                                                      index=0,
                                                                      flattened_idx_offset=flattened_idx_offset,
                                                                      encoding_nn_output_size=self.encoding_nn_output_size)

            zero = tf.constant(0)
            input_begin = tf.pack([zero, step, zero])
            step_feature = tf.squeeze(tf.slice(feature_pl, input_begin, [-1, 1, -1]), squeeze_dims=[1])

            inputs = tf.concat(1, [step_contextual_features, step_feature])
            encodings = Network.apply_EncodingNN(inputs, FLAGS.activation_type)

            molecule_encoding = tf.expand_dims(tf.reduce_sum(encodings, 0), 0)

        with tf.variable_scope("OutputNN", reuse=True) as scope:

            prediction_op = Network.apply_OutputNN(molecule_encoding,FLAGS.activation_type)

        return prediction_op

    @staticmethod
    def get_activation_fun(activation_fun):
        '''
        :param activation_fun:
        :return:
        '''
        # logging.info("Activation fun {:}".format(activation_fun))
        if activation_fun == 'tanh':
            return tf.nn.tanh
        elif activation_fun == 'relu6':
            return tf.nn.relu6
        elif activation_fun == 'crelu':
            return tf.nn.crelu
        elif activation_fun == 'relu':
            return tf.nn.relu
        elif activation_fun == 'identity':
            logging.info("Activation fun {:}".format(activation_fun))
            return tf.identity
        else:
            raise Exception("Inavlid activation function {}".format(activation_fun))

    @staticmethod
    def get_initializer(initializer):
        '''
        :param initializer:
        :return:
        '''
        logging.info("Weight Initializer:{}".format(initializer))
        if initializer == 'xavier':
            return tf.contrib.layers.xavier_initializer(uniform=False)
        elif initializer =='random':
            return tf.random_normal_initializer()
        elif initializer == 'one':
            return tf.ones_initializer()
        else:
            raise Exception("Inavlid initializer{}".format(initializer))

    @staticmethod
    def cond(sequence_len,
             step,
             feature_pl,
             path_pl,
             flattened_idx_offset,
             contextual_features,
             encoding_nn_output_size):
        return tf.less(step, sequence_len - 1)

    @staticmethod
    def body(sequence_len,
             step,
             feature_pl,
             path_pl,
             flattened_idx_offset,
             contextual_features,
             encoding_nn_output_size):
        '''
        :param sequence_len:
        :param step:
        :param feature_pl:
        :param path_pl:
        :param flattened_idx_offset:
        :param contextual_features:
        :param encoding_nn_output_size:
        :return:
        '''

        zero = tf.constant(0)
        one = tf.constant(1)
        input_begin = tf.pack([zero,step,zero])

        step_feature = tf.squeeze(tf.slice(feature_pl, input_begin, [-1, 1, -1]))

        input_idx = tf.slice(path_pl, input_begin, [-1, 1, 1])
        input_idx = tf.reshape(input_idx, [-1])

        output_begin = tf.pack([zero, step, one])
        tf.get_variable_scope().reuse_variables()

        contextual_features = tf.get_variable("contextual_features")

        step_contextual_features = Network.get_contextual_feature(contextual_features=contextual_features,
                                                                  index=input_idx,
                                                                  flattened_idx_offset=flattened_idx_offset,
                                                                  encoding_nn_output_size=encoding_nn_output_size)

        nn_inputs = tf.concat(1, [step_contextual_features, step_feature])
        updated_contextual_vectors = Network.apply_EncodingNN(nn_inputs, FLAGS.activation_type)
        output_idx = tf.squeeze(tf.slice(path_pl, output_begin, [-1, 1, 2]))

        contextual_features = Network.update_contextual_features(contextual_features=contextual_features,
                                                                 indices=output_idx,
                                                                 updates=updated_contextual_vectors,
                                                                 flattened_idx_offset=flattened_idx_offset)

        with tf.control_dependencies([contextual_features]):
            return (sequence_len,
                    step + 1,
                    feature_pl,
                    path_pl,
                    flattened_idx_offset,
                    contextual_features,
                    encoding_nn_output_size)

    def training(self):
        '''
        :return:
        '''
        def apply_gradient_clipping(gradient):
            if gradient is not None:
                return tf.mul(tf.clip_by_value(tf.abs(grad), 0.1, 1.), tf.sign(grad))
            else:
                return None

        # optimizer = tf.train.GradientDescentOptimizer(learning_rate=self.learning_rate)
        optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate, beta1=0.9, beta2=0.999, epsilon=1e-08,
                                           use_locking=False, name='Adam')

        loss_op = self.loss_op + FLAGS.weight_decay_factor*tf.add_n([tf.nn.l2_loss(v) for v in tf.get_collection('weights_decay')])

        gvs = optimizer.compute_gradients(loss_op)

        if FLAGS.clip_gradient:
            gvs = [(apply_gradient_clipping(grad), var) for grad, var in gvs]

        self.train_op = optimizer.apply_gradients(gvs)

    def loss(self):
        self.loss_op = 0

        loss_op = [tf.square(tf.sub(self.prediction_ops[i], self.target_pls[i])) for i in xrange(0, FLAGS.batch_size)]
        self.loss_op= tf.add_n(loss_op, name=None)/2
        tf.summary.scalar("loss", tf.squeeze(self.loss_op))

    @staticmethod
    def get_contextual_feature(contextual_features, index, flattened_idx_offset, encoding_nn_output_size):
        '''
        :param contextual_features:
        :param index:
        :param flattened_idx_offset:
        :param encoding_nn_output_size:
        :return:
        '''
        """
            Contextual vector is flatted array
            index is 1D index with
        """
        indices = index + flattened_idx_offset
        values = [indices, indices, indices, indices]
        indices = tf.pack(values, axis=1, name='pack')
        indices = indices + tf.constant([0, 1, 2, 3])
        indices = tf.reshape(indices, [-1])
        contextual_vector = tf.gather(contextual_features, indices)
        contextual_vector = tf.reshape(contextual_vector, [-1, 4 * encoding_nn_output_size])
        return contextual_vector

    @staticmethod
    def update_contextual_features(contextual_features, indices, updates, flattened_idx_offset):
        '''
        :param contextual_features:
        :param indices:
        :param updates:
        :param flattened_idx_offset:
        :return:
        '''
        first_indices, second_indices = tf.split(1, 2, indices)
        indices = tf.squeeze(first_indices + second_indices)
        indices = indices + flattened_idx_offset
        contextual_features = tf.scatter_add(contextual_features, indices, updates, use_locking=None)
        return contextual_features

    @staticmethod
    def weight_variable(shape, initializer, collection=None):

        weights =  tf.get_variable(name="weights",
                               shape=shape,
                               initializer=initializer,
                               trainable=True,
                               collections=['variables', collection])
        # Network.variable_summaries(weights, "weights")
        return weights

    @staticmethod
    def bias_variable(shape):
        """Create a bias variable with appropriate initialization."""
        biases= tf.get_variable(name="biases",
                               shape=shape,
                               initializer=tf.constant_initializer(0),
                               trainable=True)

        # Network.variable_summaries(biases, "bias")
        return biases

    @staticmethod
    def apply_EncodingNN(inputs, activation_type):
        activation_fun = Network.get_activation_fun(activation_type)
        with tf.variable_scope('hidden1') as scope:
            weights = tf.get_variable("weights")
            biases = tf.get_variable("biases")

            hidden1 = activation_fun(tf.matmul(inputs, weights) + biases)

        with tf.variable_scope('output') as scope:
            weights = tf.get_variable("weights")
            biases = tf.get_variable("biases")
            return activation_fun(tf.matmul(hidden1, weights) + biases)

    @staticmethod
    def apply_OutputNN(inputs, activation_type):
        activation_fun = Network.get_activation_fun(activation_type)
        with tf.variable_scope('hidden1') as scope:
            weights = tf.get_variable("weights")
            biases = tf.get_variable("biases")
            hidden1 = activation_fun(tf.matmul(inputs, weights) + biases)

        with tf.variable_scope('output') as scope:
            weights = tf.get_variable("weights")
            return tf.matmul(hidden1, weights)

    @staticmethod
    def variable_summaries(var, type):
        """Attach a lot of summaries to a Tensor (for TensorBoard visualization)."""
        with tf.name_scope(type):
            mean = tf.reduce_mean(var)
            tf.summary.scalar('mean', mean)
            with tf.name_scope('stddev'):
                stddev = tf.sqrt(tf.reduce_mean(tf.square(var - mean)))
            tf.summary.scalar('stddev', stddev)
            tf.summary.scalar('max', tf.reduce_max(var))
            tf.summary.scalar('min', tf.reduce_min(var))
            tf.summary.histogram('histogram', var)
