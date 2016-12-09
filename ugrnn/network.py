from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import logging
import tensorflow as tf
from ugrnn.molecule import Molecule

logger = logging.getLogger(__name__)

FLAGS = None


class Network(object):
    def __init__(self, name, encoding_nn_hidden_size, encoding_nn_output_size,
                 output_nn_hidden_size, max_seq_len, feature_pl, path_pl, sequence_len,
                 target_pl, learning_rate, initializer):
        """Build the ugrnn model up to where it may be used for inference."""
        self.name = name
        self.encoding_nn_output_size = encoding_nn_output_size
        self.output_nn_hidden_size = output_nn_hidden_size
        self.feature_pl = feature_pl
        self.path_pl = path_pl
        self.target_pl = target_pl
        self.sequence_len = sequence_len
        self.max_seq_len = max_seq_len
        self.encoding_nn_hidden_size = encoding_nn_hidden_size
        self.flattened_idx_offset = tf.range(0, sequence_len) * max_seq_len * 4
        self.encoding_nn_input_size = 4 * encoding_nn_output_size + Molecule.num_of_features()
        self.learning_rate = learning_rate
        self.initializer_fun = Network.get_initializer(initializer)

        self.inference()
        self.loss()
        self.training(FLAGS.clip_gradient)

    def inference(self):
        '''
        :return:
        '''
        with tf.variable_scope("EncodingNN") as scope:
            step = tf.constant(0)
            contextual_features = tf.get_variable("contextual_features",
                                                  [self.max_seq_len * self.max_seq_len * 4,
                                                   self.encoding_nn_output_size],
                                                  dtype=tf.float32,
                                                  initializer=tf.constant_initializer(0),
                                                  trainable=False)

            contextual_features = contextual_features.assign(
                tf.zeros([self.max_seq_len * self.max_seq_len * 4, self.encoding_nn_output_size],
                         dtype=tf.float32))

            with tf.variable_scope('hidden1') as scope:
                weights = Network.weight_variable([self.encoding_nn_input_size, self.encoding_nn_hidden_size],
                                                  initializer=self.initializer_fun)
                Network.variable_summaries(weights, scope.name + '/weights')

                biases = Network.bias_variable([self.encoding_nn_hidden_size])
                Network.variable_summaries(biases, scope.name + '/biases')

            with tf.variable_scope('output') as scope:
                weights = Network.weight_variable([self.encoding_nn_hidden_size, self.encoding_nn_output_size],
                                                  initializer=self.initializer_fun)
                Network.variable_summaries(weights, scope.name + '/weights')

                biases = Network.bias_variable([self.encoding_nn_output_size])
                Network.variable_summaries(biases, scope.name + '/biases')

            _, step, _, _, _, contextual_features, _ = tf.while_loop(Network.cond, Network.body,
                                                                     [self.sequence_len, step, self.feature_pl,
                                                                      self.path_pl,
                                                                      self.flattened_idx_offset, contextual_features,
                                                                      self.encoding_nn_output_size],
                                                                     back_prop=True,
                                                                     swap_memory=False, name=None)

            # use flattened indices1
            step_contextual_features = Network.get_contextual_feature(contextual_features=contextual_features,
                                                                      index=0,
                                                                      flattened_idx_offset=self.flattened_idx_offset,
                                                                      encoding_nn_output_size=self.encoding_nn_output_size)

            input_begin = tf.get_variable("input_begin", [3], dtype=tf.int32)
            input_begin = tf.scatter_update(input_begin, 1, step, use_locking=None)
            step_feature = tf.squeeze(tf.slice(self.feature_pl, input_begin, [-1, 1, -1]), squeeze_dims=[1])

            inputs = tf.concat(1, [step_contextual_features, step_feature])
            encodings = Network.apply_EncodingNN(inputs, FLAGS.activation_type)

            molecule_encoding = tf.expand_dims(tf.reduce_sum(encodings, 0), 0)

        with tf.variable_scope("OutputNN") as scope:
            hidden1 = Network.nn_layer(input_tensor=molecule_encoding,
                                       input_dim=self.encoding_nn_output_size,
                                       output_dim=self.output_nn_hidden_size,
                                       layer_name='hidden1',
                                       act=Network.get_activation_fun(FLAGS.activation_type),
                                       initializer=self.initializer_fun)

            self.prediction_op = Network.nn_layer(input_tensor=hidden1,
                                                  input_dim=self.output_nn_hidden_size,
                                                  output_dim=1,
                                                  layer_name='output',
                                                  act=tf.identity,
                                                  initializer=self.initializer_fun)

    @staticmethod
    def get_activation_fun(activation_fun):
        '''
        :param activation_fun:
        :return:
        '''
        logging.info("Activation fun {:}".format(activation_fun))
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
        if initializer == 'xavier':
            logging.info("Weight Initializer:{}".format(initializer))
            return tf.contrib.layers.xavier_initializer
        elif initializer == 'one':
            logging.info("Weight Initializer:{}".format(initializer))
            return tf.ones_initializer
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
        input_begin = tf.get_variable("input_begin", [3], dtype=tf.int32, initializer=tf.constant_initializer(0),
                                      trainable=False)
        input_begin = tf.scatter_update(input_begin, 1, step, use_locking=None)

        step_feature = tf.squeeze(tf.slice(feature_pl, input_begin, [-1, 1, -1]))

        input_idx = tf.slice(path_pl, input_begin, [-1, 1, 1])
        input_idx = tf.reshape(input_idx, [-1])

        output_begin = tf.get_variable("ouput_begin", [3], dtype=tf.int32, initializer=tf.constant_initializer(0),
                                       trainable=False)
        output_begin = tf.scatter_update(output_begin, 1, step, use_locking=None)
        output_begin = tf.scatter_update(output_begin, 2, 1, use_locking=None)

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

    def training(self, clip_gradient=True):
        '''
        :return:
        '''
        def apply_gradient_clipping(gradient):
            if gradient is not None:
                return tf.mul(tf.clip_by_value(tf.abs(grad), 0.1, 1.), tf.sign(grad))
            else:
                return None

        optimizer = tf.train.GradientDescentOptimizer(learning_rate=self.learning_rate)
        # optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate, beta1=0.9, beta2=0.999, epsilon=1e-08,
        #                                    use_locking=False, name='Adam')
        if clip_gradient:
            gvs = optimizer.compute_gradients(self.loss_op)
            capped_gvs = [(apply_gradient_clipping(grad), var) for grad, var in gvs]
            self.train_op = optimizer.apply_gradients(capped_gvs)
        else:
            self.train_op = optimizer.minimize(self.loss_op)

    def loss(self):
        self.loss_op = tf.sqrt(tf.reduce_mean(tf.square(tf.sub(self.prediction_op, self.target_pl))))/2

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
    def nn_layer(input_tensor, input_dim, output_dim, layer_name, act, initializer):
        '''
        :param input_tensor:
        :param input_dim:
        :param output_dim:
        :param layer_name:
        :param act:
        :param initializer:
        :return:
        '''
        """Reusable code for making a simple neural net layer.
        It does a matrix multiply, bias add, and then uses relu to nonlinearize.
        It also sets up name scoping so that the resultant graph is easy to read,
        and adds a number of summary ops.
        """
        # Adding a name scope ensures logical grouping of the layers in the graph.
        with tf.variable_scope(layer_name) as scope:
            weights = Network.weight_variable([input_dim, output_dim], initializer)
            Network.variable_summaries(weights, layer_name + '/weights')

            biases = Network.bias_variable([output_dim])
            Network.variable_summaries(biases, layer_name + '/biases')

            preactivate = tf.matmul(input_tensor, weights) + biases
            tf.histogram_summary(layer_name + '/pre_activations', preactivate)

            activations = act(preactivate, name='activation')
            tf.histogram_summary(layer_name + '/activations', activations)

            return activations

    @staticmethod
    def weight_variable(shape, initializer):
        return tf.get_variable(name="weights",
                               shape=shape,
                               initializer=initializer(),
                               trainable=True)

    @staticmethod
    def bias_variable(shape):
        """Create a bias variable with appropriate initialization."""
        return tf.get_variable(name="biases",
                               shape=shape,
                               initializer=tf.constant_initializer(0),
                               trainable=True)

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
    def variable_summaries(var, name):
        """Attach a lot of summaries to a Tensor."""
        with tf.name_scope('summaries'):
            mean = tf.reduce_mean(var)
            tf.scalar_summary('mean/' + name, mean)
            with tf.name_scope('stddev'):
                stddev = tf.sqrt(tf.reduce_mean(tf.square(var - mean)))
            tf.scalar_summary('stddev/' + name, stddev)
            tf.scalar_summary('max/' + name, tf.reduce_max(var))
            tf.scalar_summary('min/' + name, tf.reduce_min(var))
            tf.histogram_summary(name, var)
