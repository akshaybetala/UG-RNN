from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import logging
import os

import numpy as np

from ugrnn.input_data import DataSet
from ugrnn.ugrnn import UGRNN
from ugrnn.utils import model_params
np.set_printoptions(threshold=np.inf, precision=4)

import tensorflow as tf

FLAGS = None


def main(_):
    model_dir = os.path.join(FLAGS.output_dir, FLAGS.model_name)

    if tf.gfile.Exists(model_dir):
        tf.gfile.DeleteRecursively(model_dir)
    tf.gfile.MakeDirs(model_dir)

    with tf.Graph().as_default():
        # Create a session for running Ops on the Graph.
        sess = tf.Session()

        logp_col_name = FLAGS.logp_col if FLAGS.add_logp else None

        logger.info('Loading Training dataset from {:}'.format(FLAGS.training_file))
        train_dataset = DataSet(csv_file_path=FLAGS.training_file, smile_col_name=FLAGS.smile_col,
                                target_col_name=FLAGS.target_col, logp_col_name=logp_col_name,
                                contract_rings=FLAGS.contract_rings)

        logger.info('Loading validation dataset from {:}'.format(FLAGS.validation_file))
        validation_dataset = DataSet(csv_file_path=FLAGS.validation_file, smile_col_name=FLAGS.smile_col,
                                     target_col_name=FLAGS.target_col, logp_col_name=logp_col_name,
                                     contract_rings=FLAGS.contract_rings)

        logger.info("Creating Graph.")

        ugrnn_model = UGRNN(FLAGS.model_name, encoding_nn_hidden_size=FLAGS.model_params[0],
                            encoding_nn_output_size=FLAGS.model_params[1], output_nn_hidden_size=FLAGS.model_params[2],
                            batch_size=FLAGS.batch_size, learning_rate=0.001, add_logp=FLAGS.add_logp)

        logger.info("Succesfully created graph.")

        init = tf.global_variables_initializer()
        sess.run(init)
        logger.info('Run the Op to initialize the variables')
        ugrnn_model.train(sess, FLAGS.max_epochs, train_dataset, validation_dataset, model_dir)
        ugrnn_model.save_model(sess, model_dir, FLAGS.max_epochs)


if __name__ == '__main__':
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser()

    parser.add_argument('--model_name', type=str, default='default_model',
                        help='Name of the model')

    parser.add_argument('--max_epochs', type=int, default=200,
                        help='Number of epochs to run trainer.')

    parser.add_argument('--batch_size', type=int, default=10,
                        help='Batch size.')

    parser.add_argument('--model_params', help="Model Parameters", dest="model_params", type=model_params)

    parser.add_argument('--learning_rate', type=float, default=0.001,
                        help='Initial learning rate')

    parser.add_argument('--output_dir', type=str, default='train',
                        help='Directory for storing the trained models')

    parser.add_argument('--training_file', type=str, default='ugrnn/data/delaney/train_delaney.csv',
                        help='Path to the csv file containing training data set')

    parser.add_argument('--validation_file', type=str, default='ugrnn/data/delaney/validate_delaney.csv',
                        help='Path to the csv file containing validation data set')

    parser.add_argument('--smile_col', type=str, default='smiles')

    parser.add_argument('--logp_col', type=str, default='logp')

    parser.add_argument('--target_col', type=str, default='solubility')

    parser.add_argument('--contract_rings', dest='contract_rings',
                        action='store_true')
    parser.set_defaults(contract_rings=False)

    parser.add_argument('--add_logp', dest='add_logp',
                        action='store_true')
    parser.set_defaults(add_logp=False)

    parser.add_argument('--clip_gradient', dest='clip_gradient',
                        action='store_true')
    parser.set_defaults(clip_gradient=False)

    FLAGS = parser.parse_args()

    tf.app.run(main=main)
