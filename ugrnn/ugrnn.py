"""Trains and Evaluates the model."""
# pylint: disable=missing-docstring
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math

import numpy as np
import logging
from ugrnn import network, input_data
from ugrnn.molecule import Molecule
from ugrnn.loss import Loss

np.set_printoptions(threshold=np.inf)

from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf

import argparse
import matplotlib.pyplot as plt

# Basic model parameters as external flags.
FLAGS = None


class UGRNN(object):
    model_types = [(7, 3, 5),
                   (7, 4, 5),
                   (7, 5, 5),
                   (7, 6, 5),
                   (7, 7, 5),
                   (7, 8, 5),
                   (7, 9, 5),
                   (7, 10, 5),
                   (7, 11, 5),
                   (7, 12, 5),
                   (3, 3, 5),
                   (4, 3, 5),
                   (5, 3, 5),
                   (6, 3, 5),
                   (7, 3, 5),
                   (8, 3, 5),
                   (9, 3, 5),
                   (10, 3, 5),
                   (11, 3, 5),
                   (12, 3, 5)]

    def __init__(self, sess):
        logger.info("Creating the Network")

        self.feature_pl = tf.placeholder(tf.float32, shape=[None, None, Molecule.num_of_features()])
        self.path_pl = tf.placeholder(tf.int32, shape=[None, None, 3])
        self.target_pl = tf.placeholder(tf.float32)
        self.sequence_len_pl = tf.placeholder(tf.int32)
        self.global_step = tf.Variable(0, name='global_step', trainable=False)
        self.global_step_update_op = tf.assign(self.global_step, tf.add(self.global_step, tf.constant(1)))
        self.no_of_models = 1

        self.learning_rate = tf.train.exponential_decay(learning_rate=FLAGS.learning_rate,
                                                        global_step=self.global_step,
                                                        decay_steps=100,
                                                        decay_rate=.96,
                                                        staircase=False)

        # self.learning_rate = tf.maximum(self.learning_rate, FLAGS.min_learning_rate)

        logger.info('Initial learning rate: {:}'.format(FLAGS.learning_rate))

        self.prediction_ops = []
        self.loss_ops = []
        self.train_ops = []
        self.sess = sess

        logger.info("Ensemble {:}".format(FLAGS.ensemble))

        if FLAGS.ensemble:

            logger.info("Total No of Models {:}".format(FLAGS.ensemble))
            self.no_of_models = FLAGS.ensemble_size

            for i in xrange(self.no_of_models):
                model_name = "model_{}".format(i)
                logger.info("Building model {:} with parameters {:}".format(model_name, self.model_types[i]))
                encoding_nn_hidden_size = self.model_types[i][0]
                encoding_nn_output_size = self.model_types[i][1]
                output_nn_hidden_size = self.model_types[i][2]
                self.add_model(model_name, encoding_nn_hidden_size, encoding_nn_output_size, output_nn_hidden_size)

        else:
            model_name = "model_{}".format(FLAGS.model_no)
            logger.info("Building model {:} with parameters {:}".format(model_name, self.model_types[FLAGS.model_no]))
            encoding_nn_hidden_size = self.model_types[FLAGS.model_no][0]
            encoding_nn_output_size = self.model_types[FLAGS.model_no][1]
            output_nn_hidden_size = self.model_types[FLAGS.model_no][2]
            self.add_model(model_name, encoding_nn_hidden_size, encoding_nn_output_size, output_nn_hidden_size)

    def add_model(self, model_name, encoding_nn_hidden_size, encoding_nn_output_size, output_nn_hidden_size):

        with tf.variable_scope(model_name) as scope:
            # Build a Graph that computes predictions from the inference model.
            model = network.Network(name=model_name,
                                    encoding_nn_hidden_size=encoding_nn_hidden_size,
                                    encoding_nn_output_size=encoding_nn_output_size,
                                    output_nn_hidden_size=output_nn_hidden_size,
                                    feature_pl=self.feature_pl,
                                    path_pl=self.path_pl,
                                    target_pl=self.target_pl,
                                    sequence_len=self.sequence_len_pl,
                                    learning_rate=self.learning_rate,
                                    loss_type=FLAGS.loss_type,
                                    activation_type=FLAGS.activation_fun,
                                    max_seq_len=FLAGS.max_seq_len)

            self.prediction_ops.append(model.prediction_op)
            self.loss_ops.append(model.loss_op)
            self.train_ops.append(model.train_op)

    def optimize(self, dataset):
        n = UGRNN.num_of_networks
        total_loss_values = np.zeros(n)
        while dataset.epochs_completed < 1:
            feed_dict = self.fill_feed_dict(dataset)
            loss_values, prediction_values = self.sess.run([self.loss_ops, self.prediction_ops],
                                                           feed_dict=feed_dict)
            total_loss_values += np.array(loss_values)

        # Get the 10 netowrks with minimum error
        index_of_best_networks = total_loss_values.argsort()[int(-n / 2):]
        self.prediction_ops = [self.prediction_ops[index] for index in index_of_best_networks]

    def train(self, dataset, epochs=1):
        for epoch in xrange(epochs):
            dataset.reset_epoch()
            for i in xrange(dataset.num_examples):
                feed_dict = self.fill_feed_dict(dataset)
                _ = self.sess.run([self.train_ops], feed_dict=feed_dict)
            self.sess.run([self.global_step_update_op])

    def loss(self, dataset):
        predictions = self.predict(dataset)
        tragets = dataset.labels
        error = Loss.get_error(FLAGS.loss_type, predictions, tragets)
        return error

    def get_learning_rate(self):
        return self.sess.run([self.learning_rate])

    def predict(self, dataset):
        dataset.reset_epoch()
        predictions = []
        i=0
        while dataset.epochs_completed < 1:
            feed_dict = self.fill_feed_dict(dataset)
            prediction_values = self.sess.run([self.prediction_ops], feed_dict=feed_dict)
            predictions.append(np.mean(prediction_values))
            # print(i, prediction_values)
            i+=1

        return np.array(predictions)

    def fill_feed_dict(self, dataset):
        molecules_feed, targets_feed = dataset.next_molecule()

        feed_dict = {self.feature_pl: molecules_feed.feature_vector,
                     self.path_pl: molecules_feed.directed_graphs,
                     self.target_pl: targets_feed,
                     self.sequence_len_pl: molecules_feed.feature_vector.shape[1]}

        return feed_dict


def main(_):
    with tf.Graph().as_default():
        # Create a session for running Ops on the Graph.
        sess = tf.Session()

        logger.info("Creating Graph.")
        ugrnn_model = UGRNN(sess=sess)
        logger.info("Succesfully created graph.")

        logger.info('Run the Op to initialize the variables')
        init = tf.initialize_all_variables()
        sess.run(init)

        logger.info('Loading {:} dataset'.format(FLAGS.dataset))
        logger.info('Contract Rings: {:}'.format(FLAGS.contract_rings))

        data_sets = input_data.read_data_sets(FLAGS.dataset, FLAGS.contract_rings )
        train_dataset = data_sets.train
        validation_dataset = data_sets.validation
        test_dataset = data_sets.test

        logger.info('Start Training')

        for epoch in xrange(0, FLAGS.max_epochs):
            ugrnn_model.train(train_dataset)
            if epoch % 1 == 0:
                train_error = ugrnn_model.loss(train_dataset)
                validation_error = ugrnn_model.loss(validation_dataset)
                plt.scatter(epoch, train_error)
                learning_rate = ugrnn_model.get_learning_rate()
                plt.pause(0.05)
                logger.info(
                    "Epoch: {:}, Learning rate {:} Train Loss: {:}, Validation Loss {:}".format(epoch, learning_rate,
                                                                                                train_error,
                                                                                                validation_error))

        logger.info('Training Finished')

        logger.info('Optimize network')
        ugrnn_model.optimize(data_sets.validation)
        predictions = ugrnn_model.predict(test_dataset)
        error = Loss.get_error(FLAGS.loss_type, test_dataset.labels, predictions)
        logger.info("Loss: {:}".format(error))


if __name__ == '__main__':
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)
    logger = logging.getLogger(__name__)

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

    parser.add_argument('--dataset', type=str, default='delaney',
                        help='The data used for testing and training')

    parser.add_argument('--summaries_dir', type=str, default='/tmp/ugrnn_logs',
                        help='Summaries directory')

    parser.add_argument('--ensemble', type=bool, default=False,
                        help='The model used from the ensemble')

    parser.add_argument('--ensemble_size', type=int, default=10,
                        help='Count of models used in the ensemble')

    parser.add_argument('--model_no', type=int, default=0,
                        help='The model_type you want used to build the netwrok')

    parser.add_argument('--save_result', type=bool, default=False,
                        help='Save result')

    parser.add_argument('--contract_rings', dest='contract_rings', action='store_true')
    parser.add_argument('--no-contract_rings', dest='contract_rings', action='store_false')
    parser.set_defaults(contract_rings=False)

    FLAGS = parser.parse_args()
    plt.axis([0, FLAGS.max_epochs, 0, 4])
    plt.ion()
    tf.app.run()
