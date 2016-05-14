"""Trains and Evaluates the model."""
# pylint: disable=missing-docstring
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time
import input_data
import ugrnn
import utils


from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf


# Basic model parameters as external flags.
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_float('learning_rate', 0.01, 'Initial learning rate.')
flags.DEFINE_integer('max_steps', 80000, 'Number of steps to run trainer.')
flags.DEFINE_integer('batch_size',1, 'Number of steps to run trainer.')
flags.DEFINE_string('train_dir', 'train', 'Directory to put the training data.')
flags.DEFINE_integer('initial_feature_vector_size',utils.num_of_features(),'Size of the individual feature for all the nodes' )
flags.DEFINE_integer('max_seq_len',100,'Size of the maximum molecule')

def placeholder_inputs():
  """Generate pl variables to represent the input tensors.

  These pls are used as inputs by the rest of the model building
  code and will be fed from the downloaded data in the .run() loop, below.

  Args:
    batch_size: The batch size will be baked into both pls.

  Returns:
    molecule_pl: Images pl.
    targets_pl: targets pl.
  """
  # Note that the shapes of the pls match the shapes of the full
  # image and target tensors, except the first dimension is now batch_size
  # rather than the full size of the train or test data sets.

  feature_pl = tf.placeholder(tf.float32, shape=[None,None,FLAGS.initial_feature_vector_size])
  path_pl = tf.placeholder(tf.int32, shape=[None,None,2])
  targets_pl = tf.placeholder(tf.float32)
  sequence_len_pl = tf.placeholder(tf.int32)
  return feature_pl,path_pl, targets_pl, sequence_len_pl


def fill_feed_dict(data_set, 
                  feature_pl, 
                  path_pl, 
                  targets_pl,
                  sequence_len_pl):
  """Fills the feed_dict for training the given step.

  A feed_dict takes the form of:
  feed_dict = {
      <pl>: <tensor of values to be passed for pl>,
      ....
  }

  Args:
    data_set: The set of images and targets, from input_data.read_data_sets()
    images_pl: The images pl, from pl_inputs().
    targets_pl: The targets pl, from pl_inputs().

  Returns:
    feed_dict: The feed dictionary mapping from pls to values.
  """
  # Create the feed_dict for the pls filled with the next
  # `batch size ` examples.
  molecules_feed, targets_feed = data_set.next_batch(1)
  # print(molecules_feed[0].feature_vector.shape)
  feed_dict  = {feature_pl : molecules_feed[0].feature_vector,
                path_pl : molecules_feed[0].directed_graphs,
                targets_pl : targets_feed[0],
                sequence_len_pl : molecules_feed[0].feature_vector.shape[1]}
                
  return feed_dict


def do_eval(sess,
            loss,
            feature_pl, 
            path_pl,
            targets_pl,
            sequence_len_pl,
            data_set):
  """Runs one evaluation against the full epoch of data.

  Args:
    sess: The session in which the model has been trained.
    eval_correct: The Tensor that returns the number of correct predictions.
    images_pl: The images pl.
    targets_pl: The targets pl.
    data_set: The set of images and targets to evaluate, from
      input_data.read_data_sets().
  """
  # And run one epoch of eval.
  true_loss = 0  # Counts the number of correct predictions.
  steps_per_epoch = data_set.num_examples // FLAGS.batch_size
  num_examples = steps_per_epoch * FLAGS.batch_size
  for step in xrange(steps_per_epoch):
    feed_dict = fill_feed_dict(data_set,
                               feature_pl, 
                               path_pl,
                               targets_pl,
                               sequence_len_pl)
    true_loss += sess.run(loss, feed_dict=feed_dict)
    # print(true_loss)
  precision = true_loss / num_examples
  print('  Num examples: %d Loss @ 1: %0.04f' %
        (num_examples, precision))


def do_prediction(sess,
            predictions,
            feature_pl, 
            path_pl,
            targets_pl,
            sequence_len_pl,
            data_set):
  """Runs one evaluation against the full epoch of data.

  Args:
    sess: The session in which the model has been trained.
    eval_correct: The Tensor that returns the number of correct predictions.
    images_pl: The images pl.
    targets_pl: The targets pl.
    data_set: The set of images and targets to evaluate, from
      input_data.read_data_sets().
  """
  # And run one epoch of eval.
  true_loss = 0  # Counts the number of correct predictions.
  steps_per_epoch = data_set.num_examples // FLAGS.batch_size
  num_examples = steps_per_epoch * FLAGS.batch_size
  result = []
  for step in xrange(steps_per_epoch):
    feed_dict = fill_feed_dict(data_set,
                               feature_pl, 
                               path_pl,
                               targets_pl,
                               sequence_len_pl)
    sess.run(predictions, feed_dict=feed_dict)


def run_training(nn1_hidden_size,nn1_output_size,nn2_hidden_size):
  """Train ugrnn for a number of steps."""
  # Get the sets of molecules and targets for training, validation, and
  # test.
  data_sets = input_data.read_data_sets()
  
  # Tell TensorFlow that the model will be built into the default Graph.
  with tf.Graph().as_default():

    # Create a session for running Ops on the Graph.
    sess = tf.Session()

    # Generate pls for the images and targets.
    feature_pl,path_pl, targets_pl,sequence_len_pl = placeholder_inputs()

    predictions = []

    for i in xrange(len(nn1_hidden_size)):
      with tf.variable_scope("model{}".format(i)) as scope:
        # Build a Graph that computes predictions from the inference model.
        prediction = ugrnn.create_model(nn1_hidden_size[i],
                                          nn1_output_size[i],
                                          nn2_hidden_size[i],
                                          feature_pl,
                                          path_pl,
                                          sequence_len_pl)
        predictions.append(prediction)

    prediction = tf.add_n(predictions)/len(nn1_hidden_size)

    # # Add to the Graph the Ops for loss calculation.
    loss = ugrnn.loss(prediction, targets_pl)

    # # # # Add to the Graph the Ops that calculate and apply gradients.
    train_op = ugrnn.training(loss, FLAGS.learning_rate)

    # Build the summary operation based on the TF collection of Summaries.
    summary_op = tf.merge_all_summaries()

    # Create a saver for writing training checkpoints.
    saver = tf.train.Saver()
    
    # Instantiate a SummaryWriter to output summaries and the Graph.
    summary_writer = tf.train.SummaryWriter(FLAGS.train_dir, sess.graph)

    # Run the Op to initialize the variables.
    init = tf.initialize_all_variables()
    sess.run(init)


    # And then after everything is built, start the training loop.
    for step in xrange(FLAGS.max_steps):
      start_time = time.time()

      # Fill a feed dictionary with the actual set of images and targets
      # for this particular training step.
      feed_dict = fill_feed_dict(data_sets.train,
                                 feature_pl, 
                                 path_pl,
                                 targets_pl,
                                 sequence_len_pl)

      
      # Run one step of the model.  The return values are the activations
      # from the `train_op` (which is discarded) and the `loss` Op.  To
      # inspect the values of your Ops or variables, you may include them
      # in the list passed to sess.run() and the value tensors will be
      # returned in the tuple from the call.

      _,loss_value,predict = sess.run([train_op,loss,prediction],
                               feed_dict=feed_dict)
      
      duration = time.time() - start_time
      # Write the summaries and print an overview fairly often.
      # if step % 100 == 0:
      # #   # Print status to stdout.
      #   # print('Step %d: loss = %.2f (%.3f sec)' % (step, loss_value, duration))
      # #   # Update the events file.
      #   summary_str = sess.run(summary_op, feed_dict=feed_dict)
      #   summary_writer.add_summary(summary_str, step)
      #   summary_writer.flush()

      # Save a checkpoint and evaluate the model periodically.
      if (step + 1) % 3200 == 0 or (step + 1) == FLAGS.max_steps:
        saver.save(sess, FLAGS.train_dir, global_step=step)
        # Evaluate against the training set.
        print('Training Data Eval:')
        do_eval(sess,
                loss,
                feature_pl, 
                path_pl,
                targets_pl,
                sequence_len_pl,
                data_sets.train)
        # Evaluate against the validation set.
        print('Validation Data Eval:')
        do_eval(sess,
                loss,
                feature_pl, 
                path_pl,
                targets_pl,
                sequence_len_pl,
                data_sets.validation)
        # Evaluate against the test set.
        print('Test Data Eval:')
        do_eval(sess,
                loss,
                feature_pl, 
                path_pl,
                targets_pl,
                sequence_len_pl,
                data_sets.test)
        print("=====================================================================")

    do_prediction(sess,
                prediction,
                feature_pl, 
                path_pl,
                targets_pl,
                sequence_len_pl,
                data_sets.test)

def main(_):
  nn1_hidden_size = [7,7,7,7,7,7,7,7,7,7,3,4,5,6,7,8,9,10,11,12] 
  nn1_output_size = [3,4,5,6,7,8,9,10,11,12,3,3,3,3,3,3,3,3,3,3]
  nn2_hidden_size = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
  
  run_training(nn1_hidden_size,nn1_output_size,nn2_hidden_size)


if __name__ == '__main__':
  tf.app.run()
