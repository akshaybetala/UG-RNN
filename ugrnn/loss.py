import numpy as np
import tensorflow as tf


class Loss(object):
    @staticmethod
    def get_loss_ops(loss_type, predictions, targets):
        """

        :param loss_type:
        :param predictions:
        :param targets:
        :return:
        """
        if loss_type == 'rmse':
            return Loss.rmse_loss_ops(predictions, targets)
        elif loss_type == 'aae':
            return Loss.aae_loss_ops(predictions, targets)
        else:
            raise Exception("Invalid Loss type {}".format(loss_type))

    @staticmethod
    def rmse_loss_ops(predictions, targets):
        return tf.nn.l2_loss(tf.sub(predictions, targets)) / 2

    @staticmethod
    def aae_loss_ops(predictions, targets):
        return tf.reduce_sum(tf.abs(predictions - targets), name='l1_loss') / 2

    @staticmethod
    def get_error(loss_type, predictions, targets):
        """

        :param loss_type:
        :param predictions:
        :param targets:
        :return:
        """
        if loss_type == 'rmse':
            return Loss.rmse(predictions, targets)
        elif loss_type == 'aae':
            return Loss.aae(predictions, targets)
        else:
            raise Exception("Invalid Loss type {}".format(loss_type))

    @staticmethod
    def rmse(predictions, targets):
        return np.sqrt(((predictions - targets) ** 2).mean())

    @staticmethod
    def aae(predictions, targets):
        return np.mean(np.abs(predictions - targets))
