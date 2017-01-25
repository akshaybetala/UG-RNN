from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import logging
import os

import numpy as np

from ugrnn.input_data import DataSet
from ugrnn.ugrnn import UGRNN
from ugrnn.utils import model_params, get_metric
np.set_printoptions(threshold=np.inf)

import tensorflow as tf

FLAGS = None


def save_results(file_path, targets, predictions):
    data = np.array([targets, predictions])
    data = data.T
    f = open(file_path, 'w+')
    np.savetxt(f, data, delimiter=',', fmt=['%.4f', '%.4f'], header="Target, Prediction")
    f.close()


def get_prediction_from_model(model_name, encoding_nn_hidden_size, encoding_nn_output_size,
                              output_nn_hidden_size, test_dataset, validation_dataset):
    model_dir = os.path.join(FLAGS.output_dir, model_name)

    if not tf.gfile.Exists(model_dir):
        raise Exception("Invalid path or the model paramter doesnot exist")

    with tf.Graph().as_default():
        # Create a session for running Ops on the Graph.
        sess = tf.Session()

        logger.info("Creating Graph.")

        ugrnn_model = UGRNN(model_name, encoding_nn_hidden_size=encoding_nn_hidden_size,
                            encoding_nn_output_size=encoding_nn_output_size , output_nn_hidden_size=output_nn_hidden_size,
                            add_logp=FLAGS.add_logp)

        logger.info("Succesfully created graph.")

        init = tf.global_variables_initializer()
        sess.run(init)
        logger.info('Run the Op to initialize the variables')

        logger.info('Restoring model parameters')
        ugrnn_model.restore_model(sess, model_dir)

        prediction_validate = ugrnn_model.predict(sess, validation_dataset)
        prediction_test = ugrnn_model.predict(sess, test_dataset)

    test_results_file_path = os.path.join(model_dir, "test_result.csv")
    validation_results_file_path = os.path.join(model_dir, "validation_result.csv")

    save_results(test_results_file_path, test_dataset.labels, prediction_test)
    save_results(validation_results_file_path, validation_dataset.labels, prediction_validate)

    return prediction_validate, prediction_test


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


def main(_):
    logger.info('Loading test dataset from {:}'.format(FLAGS.test_file))
    test_dataset = DataSet(csv_file_path=FLAGS.test_file, smile_col_name=FLAGS.smile_col,
                           target_col_name=FLAGS.target_col, logp_col_name=FLAGS.logp_col,
                           contract_rings=FLAGS.contract_rings)

    logger.info('Loading validation dataset from {:}'.format(FLAGS.validation_file))
    validation_dataset = DataSet(csv_file_path=FLAGS.validation_file, smile_col_name=FLAGS.smile_col,
                                 target_col_name=FLAGS.target_col, logp_col_name=FLAGS.logp_col,
                                 contract_rings=FLAGS.contract_rings)

    for i in xrange(0,len(FLAGS.model_names)):
        prediction_validate, prediction_test = get_prediction_from_model(FLAGS.model_names[i], FLAGS.model_params[i][0],
                                                                         FLAGS.model_params[i][1], FLAGS.model_params[i][2],
                                                                         test_dataset, validation_dataset)
        prediction_metric = get_metric(prediction_validate, validation_dataset.labels)
        logger.info("RMSE: {:}, AAE: {:}".format(prediction_metric[0], prediction_metric[1]))


if __name__ == '__main__':
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser()

    parser.add_argument('--model_names', nargs='+', type=str,
                        help='Name of the models used for prediction')

    parser.add_argument('--model_params', help="Model Parameters", dest="model_params", type=model_params, nargs='+')

    parser.add_argument('--output_dir', type=str, default='train',
                        help='Root Directory where the  model parameters are stored')

    parser.add_argument('--test_file', type=str, default='ugrnn/data/delaney/validate_delaney.csv',
                        help='Path to the csv file containing test data set')

    parser.add_argument('--validation_file', type=str, default='ugrnn/data/delaney/test_delaney.csv',
                        help='Path to the csv file containing validation data set')

    parser.add_argument('--smile_col', type=str, default='smiles')

    parser.add_argument('--logp_col',  type=str, default='logp')

    parser.add_argument('--target_col', type=str, default='solubility')

    parser.add_argument('--contract_rings', dest='contract_rings',
                        action='store_true')
    parser.set_defaults(contract_rings=False)

    parser.add_argument('--add_logp', dest='add_logp',
                        action='store_true')
    parser.set_defaults(add_logp=False)

    FLAGS = parser.parse_args()
    assert len(FLAGS.model_params) == len(FLAGS.model_names)

    tf.app.run(main=main)


