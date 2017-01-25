from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from sklearn import linear_model

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
    np.savetxt(f, data, delimiter=',', fmt=['%.4f', '%.4f'], header="Target, Prediction",comments="")
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

def get_next_best_model(index, current_prediction, all_predictions, targets):
    print("======================================")
    no_of_models = len(all_predictions)

    current_error = (get_metric(current_prediction,targets))[0]
    print("Current error {:}".format(current_error))
    next_best_model_index = -1

    for i in xrange(0, no_of_models):
        temp_prediction = (index*current_prediction + all_predictions[i])/(index+1)
        metric = get_metric(temp_prediction, targets)
        print(metric[0])
        if metric[0] < current_error:
            next_best_model_index = i
            current_error = metric[0]

    return next_best_model_index

def main(_):
    logger.info('Loading test dataset from {:}'.format(FLAGS.test_file))
    test_dataset = DataSet(csv_file_path=FLAGS.test_file, smile_col_name=FLAGS.smile_col,
                           target_col_name=FLAGS.target_col, logp_col_name=FLAGS.logp_col,
                           contract_rings=FLAGS.contract_rings)

    logger.info('Loading validation dataset from {:}'.format(FLAGS.validation_file))
    validation_dataset = DataSet(csv_file_path=FLAGS.validation_file, smile_col_name=FLAGS.smile_col,
                                 target_col_name=FLAGS.target_col, logp_col_name=FLAGS.logp_col,
                                 contract_rings=FLAGS.contract_rings)

    validation_predictions = np.empty((len(FLAGS.model_names), validation_dataset.num_examples))
    test_predictions_ = np.empty((len(FLAGS.model_names), test_dataset.num_examples))

    for i in xrange(0,len(FLAGS.model_names)):
        predictions = get_prediction_from_model(FLAGS.model_names[i], FLAGS.model_params[i][0],
                                                                         FLAGS.model_params[i][1], FLAGS.model_params[i][2],
                                                                         test_dataset, validation_dataset)

        prediction_metric = get_metric(predictions[0], validation_dataset.labels)
        logger.info("Model {:} RMSE: {:}, AAE: {:}".format(FLAGS.model_names[i], prediction_metric[0], prediction_metric[1]))

        validation_predictions[i, :] = predictions[0]
        test_predictions_[i, :] = predictions[1]

    if FLAGS.optimize_ensemble:
        lr = linear_model.LinearRegression(fit_intercept=False)
        lr = lr.fit(validation_predictions.T, validation_dataset.labels)
        emsemble_preditions = lr.predict(test_predictions_.T)
    else:
        emsemble_preditions = np.mean(test_predictions_, axis=0)

    prediction_metric = get_metric(emsemble_preditions, test_dataset.labels)
    logger.info("RMSE: {:}, AAE: {:}".format(prediction_metric[0], prediction_metric[1]))

    final_prediction_path = os.path.join(FLAGS.output_dir,"ensemble_test_prediction.csv")
    save_results(final_prediction_path, test_dataset.labels, emsemble_preditions)


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

    parser.add_argument('--optimize_ensemble', dest='optimize_ensemble',
                        action='store_true')
    parser.set_defaults(optimize_ensemble=False)

    FLAGS = parser.parse_args()
    assert len(FLAGS.model_params) == len(FLAGS.model_names)

    tf.app.run(main=main)


