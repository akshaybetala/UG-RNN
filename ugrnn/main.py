"""Trains and Evaluates the model."""
# pylint: disable=missing-docstring
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import logging
import os

import numpy as np

from ugrnn import network, input_data
from ugrnn.molecule import Molecule

np.set_printoptions(threshold=np.inf)

from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf

import argparse
import matplotlib.pyplot as plt

# Basic model parameters as external flags.
FLAGS = None


class UGRNN(object):
    model_types = [(20, 20, 10),
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

        decay_rate = FLAGS.learning_rate / FLAGS.max_epochs
        decay_steps = FLAGS.max_epochs
        self.learning_rate = tf.train.exponential_decay(learning_rate=FLAGS.learning_rate,
                                                        global_step=self.global_step,
                                                        decay_steps=decay_steps,
                                                        decay_rate=decay_rate,
                                                        staircase=False)

        # self.learning_rate = tf.maximum(self.learning_rate, FLAGS.min_learning_rate)

        logger.info('Initial learning rate: {:}'.format(FLAGS.learning_rate))

        self.prediction_ops = []
        self.loss_ops = []
        self.train_ops = []
        self.sess = sess

        logger.info("Ensemble {:}".format(FLAGS.ensemble))

        if FLAGS.ensemble:
            logger.info("Total number of Models {:}".format(len(self.model_types)))
            for i, model_param in enumerate(self.model_types):
                model_name = "model_{}".format(i)
                logger.info("Building model {:} with parameters {:}".format(model_name, model_param))
                self.add_model(model_name, model_param)

        else:
            model_name = "model_{}".format(FLAGS.model_no)
            model_param = self.model_types[FLAGS.model_no]
            logger.info("Building model {:} with parameters {:}".format(model_name, model_param))
            self.add_model(model_name, model_param)

    def add_model(self, model_name, model_param):
        '''
        :param model_name:
        :return:
        '''
        encoding_nn_hidden_size = model_param[0]
        encoding_nn_output_size = model_param[1]
        output_nn_hidden_size = model_param[2]
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
                                    max_seq_len=FLAGS.max_seq_len,
                                    initializer=FLAGS.initializer)

            self.prediction_ops.append(model.prediction_op)
            self.loss_ops.append(model.loss_op)
            self.train_ops.append(model.train_op)

    def get_learning_rate(self):
        return self.sess.run([self.learning_rate])

    def optimize(self, dataset):
        '''
        :return:
        '''
        logger.info('Optimize network')
        logger.info('No of of models {:}'.format(self.no_of_models))
        self.no_of_best_models = int(self.no_of_models / 2)
        total_loss_values = np.zeros(self.no_of_models)
        while dataset.epochs_completed < 1:
            feed_dict = self.fill_feed_dict(dataset)
            loss_values, prediction_values = self.sess.run([self.loss_ops, self.prediction_ops],
                                                           feed_dict=feed_dict)
            total_loss_values += np.array(loss_values)
        index_of_best_networks = total_loss_values.argsort()[-self.no_of_best_models:]
        self.prediction_ops = [self.prediction_ops[index] for index in index_of_best_networks]

    def train(self, train_dataset, validation_dataset, epochs=1):
        '''
        :param epochs:
        :return:
        '''

        plt.subplot(2, 1, 1)
        plt.title('Training data set')
        plt.axis([0, FLAGS.max_epochs, 0, 4])

        plt.subplot(2, 1, 2)
        plt.title('Vaidation data set')
        plt.axis([0, FLAGS.max_epochs, 0, 4])

        plt.ion()
        saver = tf.train.Saver()
        logger.info('Start Training')

        for epoch in xrange(0, epochs):
            for i in xrange(train_dataset.num_examples):
                feed_dict = self.fill_feed_dict(train_dataset)
                _, _loss = self.sess.run([self.train_ops, self.loss_ops], feed_dict=feed_dict)

            train_dataset.reset_epoch(permute=True)

            self.sess.run([self.global_step_update_op])

            if epoch % 5 == 0:
                train_metric = self.evaluate(train_dataset)
                validation_metric = self.evaluate(validation_dataset, write_result=False)
                plt.subplot(2, 1, 1)
                plt.scatter(epoch, train_metric[0],  color='red', marker=".")
                plt.scatter(epoch, train_metric[1], color='blue', marker=".")

                plt.subplot(2, 1, 2)
                plt.scatter(epoch, validation_metric[0], color='red', marker=".")
                plt.scatter(epoch, validation_metric[1], color='blue', marker=".")
                learning_rate = self.get_learning_rate()
                plt.pause(0.05)
                logger.info("Epoch: {:}, Learning rate {:}  Train RMSE: {:}, Train AAE: {:} Validation RMSE {:}, Validation AAE {:}".
                            format(epoch, learning_rate[0], train_metric[0],train_metric[1], validation_metric[0], validation_metric[1]))

        if FLAGS.ensemble:
            self.optimize(validation_dataset)

        logger.info('Training Finished')

    def evaluate(self, dataset, write_result=False):
        predictions = self.predict(dataset)
        targets = dataset.labels
        if write_result:
            self.write_result(predictions, targets)
        return UGRNN.get_metric(predictions, targets)

    @staticmethod
    def get_metric(predictions, targets):
        rmse = np.sqrt(((predictions - targets) ** 2).mean())
        aae = np.mean(np.abs(predictions - targets))
        return rmse, aae

    def write_result(self, predictions, targets):
        '''
        :param predictions:
        :param targets:
        :return:
        '''
        f = open(UGRNN.get_result_file_path(FLAGS.model_no), 'w+')
        data = np.array([predictions, targets])
        data = data.T
        np.savetxt(f, data, delimiter=',', fmt=['%.4f', '%.4f'], header="Prediction, Target")
        f.close()

    def predict(self, dataset):
        '''
        :param dataset:
        :return:
        '''
        dataset.reset_epoch()
        predictions = []
        while dataset.epochs_completed < 1:
            feed_dict = self.fill_feed_dict(dataset)
            prediction_values = self.sess.run([self.prediction_ops], feed_dict=feed_dict)
            predictions.append(np.mean(prediction_values))
        return np.array(predictions)

    @staticmethod
    def get_model_file_path(model_no):
        '''
        :return:
        '''
        path = os.path.dirname(os.path.realpath(__file__))
        model_name = "model_{}".format(FLAGS.model_no)
        folder = "{:}/{:}".format(FLAGS.train_dir, model_name)
        model_file_name = "ugrnn.model" if not FLAGS.contract_rings else "ugrnn-cr.model"
        folder_path = os.path.join(path, folder)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_path = os.path.join(folder_path, model_file_name)
        return file_path

    @staticmethod
    def get_result_file_path(model_no):
        '''
        :return:
        '''
        path = os.path.dirname(os.path.realpath(__file__))
        model_name = "model_{}".format(FLAGS.model_no)
        folder = "{:}/{:}".format(FLAGS.result_dir, model_name)
        result_file_name = "ugrnn.result" if not FLAGS.contract_rings else "ugrnn-cr.result"
        folder_path = os.path.join(path, folder)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_path = os.path.join(folder_path, result_file_name)
        return file_path

    @staticmethod
    def get_log_file_path(model_no):
        '''
        :return:
        '''
        path = os.path.dirname(os.path.realpath(__file__))
        model_name = "model_{}".format(FLAGS.model_no)
        folder = "{:}/{:}".format(FLAGS.log_dir, model_name)
        log_file_name = "ugrnn.logs" if not FLAGS.contract_rings else "ugrnn-cr.logs"
        folder_path = os.path.join(path, folder)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_path = os.path.join(folder_path, log_file_name)
        return file_path

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

        logger.info('Loading {:} dataset'.format(FLAGS.dataset))
        logger.info('Contract Rings: {:}'.format(FLAGS.contract_rings))
        data_sets = input_data.read_data_sets(FLAGS.dataset, FLAGS.contract_rings)
        train_dataset = data_sets.train
        validation_dataset = data_sets.validation
        test_dataset = data_sets.test

        logger.info("Creating Graph.")
        ugrnn_model = UGRNN(sess=sess)
        logger.info("Succesfully created graph.")

        saver = tf.train.Saver()

        model_file_name = UGRNN.get_model_file_path(FLAGS.model_no)

        if FLAGS.train:
            logger.info('Run the Op to initialize the variables')
            init =tf.global_variables_initializer()
            sess.run(init)
            ugrnn_model.train(train_dataset, validation_dataset, epochs=FLAGS.max_epochs)
            saver.save(sess, model_file_name)
        else:
            init = tf.global_variables_initializer()
            sess.run(init)
            saver.restore(sess, model_file_name)
            predictions = ugrnn_model.predict(test_dataset)
            prediction_metric = ugrnn_model.evaluate(test_dataset)
            logger.info("RMSE: {:}, AAE: {:}".format(prediction_metric[0],prediction_metric[1]))

if __name__ == '__main__':
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser()

    parser.add_argument('--train', dest='train', action='store_true')
    parser.set_defaults(train=False)

    parser.add_argument('--predict', dest='predict', action='store_true')
    parser.set_defaults(predict=False)

    parser.add_argument('--max_epochs', type=int, default=2000,
                        help='Number of epochs to run trainer.')

    parser.add_argument('--learning_rate', type=float, default=0.01,
                        help='Initial learning rate')

    parser.add_argument('--learning_rate_decay_factor', type=float, default=0.98,
                        help='Initial learning rate')

    parser.add_argument('--max_seq_len', type=int, default=100,
                        help='Maximum numbers of atoms, that can be present in a Molecule')

    parser.add_argument('--train_dir', type=str, default='models',
                        help='Directory for storing the trained models')

    parser.add_argument('--activation_type', type=str, default='relu',
                        help='Summaries directory')

    parser.add_argument('--dataset', type=str, default='delaney',
                        help='The data used for testing and training')

    parser.add_argument('--initializer', type=str, default='xavier',
                        help='Initializer used to initialize the weights')

    parser.add_argument('--model_no', type=int, default=0,
                        help='Select from 0-19 to use a predefined model')

    parser.add_argument('--ensemble', dest='ensemble', action='store_true')
    parser.set_defaults(ensemble=False)

    parser.add_argument('--save_result', dest='save_result', action='store_true')
    parser.set_defaults(save_result=False)

    parser.add_argument('--contract_rings', dest='contract_rings', action='store_true')
    parser.set_defaults(contract_rings=False)

    parser.add_argument('--clip_gradient', dest='clip_gradient', action='store_true')
    parser.set_defaults(clip_gradient=False)

    FLAGS = parser.parse_args()

    if FLAGS.model_no < 0 or FLAGS.model_no > 19:
        parser.error("Invalid model number")

    if FLAGS.train and FLAGS.predict:
        parser.error('Both train and predict cannot be used together')
    if not FLAGS.train and not FLAGS.predict:
        parser.error('Either train or predict should be set')

    network.FLAGS = FLAGS
    tf.app.run()
