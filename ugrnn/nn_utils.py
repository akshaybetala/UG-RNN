import tensorflow as tf


def get_activation_fun(activation_fun):
    if activation_fun == 'tanh':
        return tf.nn.tanh
    elif activation_fun == 'relu6':
        return tf.nn.relu6
    elif activation_fun == 'crelu':
        return tf.nn.crelu
    elif activation_fun == 'relu':
        return tf.nn.relu
    elif activation_fun == 'identity':
        return tf.identity
    else:
        raise Exception(
            "Inavlid activation function {}".format(activation_fun))


def get_initializer(initializer):
    if initializer == 'xavier':
        return tf.contrib.layers.xavier_initializer(uniform=False)
    elif initializer == 'random':
        return tf.random_normal_initializer()
    elif initializer == 'one':
        return tf.ones_initializer()
    else:
        raise Exception("Inavlid initializer{}".format(initializer))


def weight_variable(shape, initializer, collection=None):
    weights = tf.get_variable(name="weights",
                              shape=shape,
                              initializer=initializer,
                              trainable=True,
                              collections=['variables', collection])
    return weights


def bias_variable(shape):
    """Create a bias variable with appropriate initialization."""
    biases = tf.get_variable(name="biases",
                             shape=shape,
                             initializer=tf.constant_initializer(0),
                             trainable=True)

    return biases
