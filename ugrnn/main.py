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
		batch_size = FLAGS.batch_size
		self.local_input_pls = \
			[tf.placeholder(tf.float32, shape=[None, None, Molecule.num_of_features()]) for i in xrange(batch_size)]
		self.path_pls = [tf.placeholder(tf.int32, shape=[None, None, 3]) for i in xrange(batch_size)]
		self.target_pls = [tf.placeholder(tf.float32) for i in xrange(batch_size)]
		self.logP_pls = [tf.placeholder(tf.float32) for i in xrange(batch_size)]
		self.sequence_len_pls = [tf.placeholder(tf.int32) for i in xrange(batch_size)]
		self.global_step = tf.Variable(0, name='global_step', trainable=False)
		self.global_step_update_op = tf.assign(self.global_step,tf.add(self.global_step, tf.constant(1)))
		self.no_of_models = 1
		self.models = []
		self.prediction_ops = []
		self.loss_ops = []
		self.train_ops = []

		self.learning_rate = FLAGS.learning_rate * tf.pow(
			FLAGS.learning_rate_decay_factor,
			tf.to_float(self.global_step), name=None)
		logger.info('Initial learning rate: {:}'.format(FLAGS.learning_rate))
		logger.info("Ensemble {:}".format(FLAGS.ensemble))

		self.sess = sess

		if FLAGS.ensemble:
			self.no_of_models = len(self.model_types)
			logger.info(
				"Total number of Models {:}".format(len(self.model_types)))
			for i, model_param in enumerate(self.model_types):
				model_name = "model_{}".format(i)
				logger.info("Building model {:} with parameters {:}".format(model_name,model_param))
				model = self.create_model(model_name, model_param)
				self.models.append(model)
				self.loss_ops.append(model.loss_op)
				self.train_ops.append(model.train_op)
				self.prediction_ops.append(model.prediction_ops[0])
		else:
			model_name = "model_{}".format(FLAGS.model_no)
			model_param = self.model_types[FLAGS.model_no]
			logger.info("Building model {:} with parameters {:}".format(model_name, model_param))
			model = self.create_model(model_name, model_param)
			self.models.append(model)
			self.loss_ops.append(model.loss_op)
			self.train_ops.append(model.train_op)
			self.prediction_ops.append(model.prediction_ops[0])
		self.summaries = tf.summary.merge_all()

	def create_model(self, model_name, model_param):
		encoding_nn_hidden_size = model_param[0]
		encoding_nn_output_size = model_param[1]
		output_nn_hidden_size = model_param[2]
		with tf.variable_scope(model_name) as scope:
			# Build a Graph that computes predictions from the inference model.
			model = network.Network(name=model_name,
									encoding_nn_hidden_size=encoding_nn_hidden_size,
									encoding_nn_output_size=encoding_nn_output_size,
									output_nn_hidden_size=output_nn_hidden_size,
									feature_pls=self.local_input_pls,
									path_pls=self.path_pls,
									logP_pls=self.logP_pls,
									target_pls=self.target_pls,
									sequence_lens=self.sequence_len_pls,
									learning_rate=self.learning_rate,
									max_seq_len=FLAGS.max_seq_len,
									initializer=FLAGS.initializer)
		return model

	def get_learning_rate(self):
		return self.sess.run([self.learning_rate])

	def save_ugrnn(self, step):
		for model in self.models:
			model.save_network(self.sess, step)

	def restore_ugrnn(self):
		for model in self.models:
			model.restore_network(self.sess)

	def optimize(self, dataset):
		logger.info('Optimize network')
		logger.info('No of of models {:}'.format(self.no_of_models))
		predictions, _ = self.predict(dataset)
		errors = []
		for i in xrange(0, self.no_of_models):
			model_prediction = predictions[:, i]
			error = self.get_metric(model_prediction, dataset.labels)
			errors.append(error[0])
		self.no_of_best_models = int(self.no_of_models / 2)
		errors = np.array(errors)
		index_of_best_networks = errors.argsort()[:self.no_of_best_models]
		logging.info(errors)
		logging.info(index_of_best_networks)
		self.prediction_ops = [self.prediction_ops[index] for index in index_of_best_networks]

	def train(self, train_dataset, validation_dataset, epochs=1):
		train_writer = tf.train.SummaryWriter(FLAGS.summaries_dir + '/train', self.sess.graph)

		plt.subplot(2, 1, 1)
		plt.title('Training data set')
		plt.axis([0, FLAGS.max_epochs, 0, 4])

		plt.subplot(2, 1, 2)
		plt.title('Vaidation data set')
		plt.axis([0, FLAGS.max_epochs, 0, 4])

		plt.ion()
		saver = tf.train.Saver()
		logger.info('Start Training')

		steps_in_epoch = train_dataset.num_examples // FLAGS.batch_size

		for epoch in xrange(0, epochs):
			for i in xrange(steps_in_epoch):
				feed_dict = self.fill_feed_dict(train_dataset,
												FLAGS.batch_size)
				_, summaries = self.sess.run([self.train_ops, self.summaries], feed_dict=feed_dict)
			train_writer.add_summary(summaries, epoch)

			train_dataset.reset_epoch(permute=True)

			self.sess.run([self.global_step_update_op])

			if epoch % 10 == 0:
				train_metric = self.evaluate(train_dataset)
				validation_metric = self.evaluate(validation_dataset)
				plt.subplot(2, 1, 1)
				plt.scatter(epoch, train_metric[0], color='red', marker=".")
				plt.scatter(epoch, train_metric[1], color='blue', marker=".")

				plt.subplot(2, 1, 2)
				plt.scatter(epoch, validation_metric[0], color='red', marker=".")
				plt.scatter(epoch, validation_metric[1], color='blue', marker=".")
				learning_rate = self.get_learning_rate()
				plt.pause(0.05)
				logger.info(
					"Epoch: {:}, Learning rate {:.8f}  Train RMSE: {:.4f}, Train AAE: {:.4f} Validation RMSE {:.4f}, Validation AAE {:.4f}".
						format(epoch, learning_rate[0], train_metric[0],
							   train_metric[1], validation_metric[0],
							   validation_metric[1],
							   precision=8))

		self.save_ugrnn(epochs)
		logger.info('Training Finished')

	def evaluate(self, dataset):
		_, predictions = self.predict(dataset)
		targets = dataset.labels
		return UGRNN.get_metric(predictions, targets)

	@staticmethod
	def get_metric(predictions, targets):
		rmse = np.sqrt(np.mean((predictions - targets) ** 2))
		aae = np.mean(np.abs(predictions - targets))
		return rmse, aae

	def predict(self, dataset):
		dataset.reset_epoch()
		individual_predictions = np.empty(
			(dataset.num_examples, len(self.prediction_ops)))
		average = []
		for i in xrange(0, dataset.num_examples):
			feed_dict = self.fill_feed_dict(dataset, 1)
			prediction_values = self.sess.run([self.prediction_ops], feed_dict=feed_dict)
			individual_predictions[i, :] = np.squeeze(prediction_values)
			average.append(np.mean(prediction_values))

		return individual_predictions, np.array(average)

	@staticmethod
	def get_result_file_path(model_no):
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
		path = os.path.dirname(os.path.realpath(__file__))
		model_name = "model_{}".format(FLAGS.model_no)
		folder = "{:}/{:}".format(FLAGS.log_dir, model_name)
		log_file_name = "ugrnn.logs" if not FLAGS.contract_rings else "ugrnn-cr.logs"
		folder_path = os.path.join(path, folder)
		if not os.path.exists(folder_path):
			os.makedirs(folder_path)
		file_path = os.path.join(folder_path, log_file_name)
		return file_path

	def fill_feed_dict(self, dataset, batch_size):
		assert batch_size <= FLAGS.batch_size
		molecules_feeds, targets_feeds = dataset.next_batch(batch_size)
		feed_dict = {}
		for i in xrange(batch_size):
			feed_dict[self.local_input_pls[i]] = molecules_feeds[i].local_input_vector
			feed_dict[self.path_pls[i]] = molecules_feeds[i].directed_graphs
			feed_dict[self.target_pls[i]] = targets_feeds[i]
			feed_dict[self.sequence_len_pls[i]] = molecules_feeds[i].local_input_vector.shape[1]

			if FLAGS.add_logP:
				feed_dict[self.logP_pls[i]] = molecules_feeds[i].logP

		return feed_dict


def main(_):
	with tf.Graph().as_default():
		# Create a session for running Ops on the Graph.
		sess = tf.Session()

		logger.info('Loading {:} dataset'.format(FLAGS.dataset))
		logger.info('Contract Rings: {:}'.format(FLAGS.contract_rings))
		data_sets = input_data.read_data_sets(FLAGS.dataset, FLAGS.add_logP, FLAGS.contract_rings)
		train_dataset = data_sets.train
		validation_dataset = data_sets.validation
		test_dataset = data_sets.test

		logger.info("Creating Graph.")
		ugrnn_model = UGRNN(sess=sess)
		logger.info("Succesfully created graph.")

		init = tf.global_variables_initializer()
		sess.run(init)

		if FLAGS.train:
			logger.info('Run the Op to initialize the variables')
			ugrnn_model.train(train_dataset, validation_dataset, epochs=FLAGS.max_epochs)
		else:
			ugrnn_model.restore_ugrnn()
			if FLAGS.ensemble:
				ugrnn_model.optimize(validation_dataset)

		_, predictions = ugrnn_model.predict(test_dataset)
		prediction_metric = ugrnn_model.get_metric(predictions, test_dataset.labels)
		logger.info("RMSE: {:}, AAE: {:}".format(prediction_metric[0], prediction_metric[1]))
		data = np.array([predictions, test_dataset.labels])
		data = data.T
		f = open("results", 'w+')
		np.savetxt(f, data, delimiter=',', fmt=['%.4f', '%.4f'], header="Prediction, Target")
		f.close()


if __name__ == '__main__':
	log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
	logging.basicConfig(level=logging.INFO, format=log_format)
	logger = logging.getLogger(__name__)

	parser = argparse.ArgumentParser()

	parser.add_argument('--train', dest='train', action='store_true')
	parser.set_defaults(train=False)

	parser.add_argument('--predict', dest='predict', action='store_true')
	parser.set_defaults(predict=False)

	parser.add_argument('--max_epochs', type=int, default=200,
						help='Number of epochs to run trainer.')

	parser.add_argument('--batch_size', type=int, default=5,
						help='Batch size.')

	parser.add_argument('--learning_rate', type=float, default=0.001,
						help='Initial learning rate')

	parser.add_argument('--learning_rate_decay_factor', type=float,
						default=0.99,
						help='Initial learning rate')

	parser.add_argument('--weight_decay_factor', type=float, default=0.0000,
						help='Weight decay factor')

	parser.add_argument('--max_seq_len', type=int, default=100,
						help='Maximum numbers of atoms, that can be present in a Molecule')

	parser.add_argument('--summaries_dir', type=str, default='train',
						help='Directory for storing graph summaries')

	parser.add_argument('--models_dir', type=str, default='models',
						help='Directory for storing the trained models')

	parser.add_argument('--activation_type', type=str, default='relu6',
						help='Summaries directory')

	parser.add_argument('--dataset', type=str, default='delaney',
						help='The data used for testing and training')

	parser.add_argument('--initializer', type=str, default='xavier',
						help='Initializer used to initialize the weights')

	parser.add_argument('--model_no', type=int, default=0,
						help='Select from 0-19 to use a predefined model')

	parser.add_argument('--ensemble', dest='ensemble', action='store_true')
	parser.set_defaults(ensemble=False)

	parser.add_argument('--save_result', dest='save_result',
						action='store_true')
	parser.set_defaults(save_result=False)

	parser.add_argument('--contract_rings', dest='contract_rings',
						action='store_true')
	parser.set_defaults(contract_rings=False)

	parser.add_argument('--add_logP', dest='add_logP',
						action='store_true')
	parser.set_defaults(add_logP=False)


	parser.add_argument('--clip_gradient', dest='clip_gradient',
						action='store_true')
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
