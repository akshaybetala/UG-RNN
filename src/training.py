"""Trains and Evaluates the model."""
# pylint: disable=missing-docstring
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time

from six.moves import xrange  # pylint: disable=redefined-builtin
import tensorflow as tf

import input_data
import ugrnn
import utils


# Basic model parameters as external flags.
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_float('learning_rate', 0.01, 'Initial learning rate.')
flags.DEFINE_integer('max_steps', 2000, 'Number of steps to run trainer.')
flags.DEFINE_integer('batch_size',1, 'Number of steps to run trainer.')
flags.DEFINE_integer('rnn_hidden_size',7, 'Number of units in hidden layer of rnn.')
flags.DEFINE_integer('rnn_output_size', 3, 'Size of the state in rnn.')
flags.DEFINE_integer('nn_hidden_size', 5, 'Number of units in hidden layer of fully connected nn.')
flags.DEFINE_string('train_dir', 'data', 'Directory to put the training data.')
flags.DEFINE_integer('initial_feature_vector_size',utils.num_of_features(),'Size of the individual feature for all the nodes' )
flags.DEFINE_integer('contextual_vector_size',FLAGS.rnn_output_size,'Size of the learned features for all nodes')
flags.DEFINE_integer('maximum_sequence_length',100,'Size of the maximum molecule')

def placeholder_inputs():
  """Generate placeholder variables to represent the input tensors.

  These placeholders are used as inputs by the rest of the model building
  code and will be fed from the downloaded data in the .run() loop, below.

  Args:
    batch_size: The batch size will be baked into both placeholders.

  Returns:
    molecule_placeholder: Images placeholder.
    targets_placeholder: targets placeholder.
  """
  # Note that the shapes of the placeholders match the shapes of the full
  # image and target tensors, except the first dimension is now batch_size
  # rather than the full size of the train or test data sets.

  feature_placeholder =[tf.placeholder(tf.int32, shape=(None,FLAGS.initial_feature_vector_size)) for i in xrange(0,FLAGS.maximum_sequence_length)]
  # feature_placeholder =tf.placeholder(tf.int32, shape=(None,None,FLAGS.initial_feature_vector_size)) 
  path_placeholder = tf.placeholder(tf.int32, shape=(None,None,2))
  targets_placeholder = tf.placeholder(tf.int32, shape=(1))
  sequence_length_placeholder =  tf.placeholder(tf.int32)
  contextual_vector_placeholder = tf.placeholder(tf.int32, shape=(None,None,FLAGS.contextual_vector_size))
  return feature_placeholder,path_placeholder, targets_placeholder, sequence_length_placeholder,contextual_vector_placeholder


def fill_feed_dict(data_set, 
                  feature_placeholder, 
                  path_placeholder, 
                  targets_placeholder, 
                  sequence_length_placeholder,
                  contextual_vector_placeholder):
  """Fills the feed_dict for training the given step.

  A feed_dict takes the form of:
  feed_dict = {
      <placeholder>: <tensor of values to be passed for placeholder>,
      ....
  }

  Args:
    data_set: The set of images and targets, from input_data.read_data_sets()
    images_pl: The images placeholder, from placeholder_inputs().
    targets_pl: The targets placeholder, from placeholder_inputs().

  Returns:
    feed_dict: The feed dictionary mapping from placeholders to values.
  """
  # Create the feed_dict for the placeholders filled with the next
  # `batch size ` examples.
  molecules_feed, targets_feed = data_set.next_batch(1)
  sequence_length = molecules_feed[0].feature_vector.shape()[0]

  feed_dict = {k: v for k, v in zip(feature_placeholder, molecules_feed[0].feature_vector)}
  feed_dict[path_placeholder] = molecules_feed[0].directed_paths
  feed_dict[targets_placeholder]= targets_feed[0]
  feed_dict[sequence_length_placeholder]=sequence_length
  feed_dict[contextual_vector_placeholder]=tf.zeros([sequence_length,sequence_length,FLAGS.contextual_vector_size])


  # feed_dict = {
  #     feature_placeholder: molecules_feed[0].feature_vector,
  #     path_placeholder: molecules_feed[0].directed_paths,
  #     targets_placeholder: targets_feed[0],
  #     sequence_length_placeholder:sequence_length,
  #     contextual_vector_placeholder:tf.zeros([sequence_length,sequence_length,FLAGS.contextual_vector_size])
  # }

  print(feed_dict)
  return feed_dict


def do_eval(sess,
            eval_correct,
            feature_placeholder, 
            path_placeholder,
            targets_placeholder,
            sequence_length_placeholder,
            contextual_vector_placeholder,
            data_set):
  """Runs one evaluation against the full epoch of data.

  Args:
    sess: The session in which the model has been trained.
    eval_correct: The Tensor that returns the number of correct predictions.
    images_placeholder: The images placeholder.
    targets_placeholder: The targets placeholder.
    data_set: The set of images and targets to evaluate, from
      input_data.read_data_sets().
  """
  # And run one epoch of eval.
  true_count = 0  # Counts the number of correct predictions.
  steps_per_epoch = data_set.num_examples // FLAGS.batch_size
  num_examples = steps_per_epoch * FLAGS.batch_size
  for step in xrange(steps_per_epoch):
    feed_dict = fill_feed_dict(data_set,
                               feature_placeholder, 
                               path_placeholder,
                               targets_placeholder,
                               sequence_length_placeholder,
                               contextual_vector_placeholder)
    true_count += sess.run(eval_correct, feed_dict=feed_dict)
  precision = true_count / num_examples
  print('  Num examples: %d  Num correct: %d  Precision @ 1: %0.04f' %
        (num_examples, true_count, precision))


def run_training():
  """Train ugrnn for a number of steps."""
  # Get the sets of molecules and targets for training, validation, and
  # test.
  data_sets = input_data.read_data_sets()

  # Tell TensorFlow that the model will be built into the default Graph.
  with tf.Graph().as_default():

    # Create a session for running Ops on the Graph.
    sess = tf.Session()

    # Generate placeholders for the images and targets.
    feature_placeholder,path_placeholder, targets_placeholder,sequence_length_placeholder,contextual_vector_placeholder = placeholder_inputs()

    # Run the Op to initialize the variables.
    init = tf.initialize_all_variables()
    sess.run(init)

    # Build a Graph that computes predictions from the inference model.
    logits = ugrnn.inference(FLAGS.rnn_hidden_size,
                            FLAGS.rnn_output_size,
                            FLAGS.nn_hidden_size,
                            feature_placeholder,
                            path_placeholder,
                            sequence_length_placeholder,
                            contextual_vector_placeholder
                            )

    # Add to the Graph the Ops for loss calculation.
    loss = ugrnn.loss(logits, targets_placeholder)

    # Add to the Graph the Ops that calculate and apply gradients.
    train_op = ugrnn.training(loss, FLAGS.learning_rate)

    # Add the Op to compare the logits to the targets during evaluation.
    eval_correct = ugrnn.evaluation(logits, targets_placeholder)

    # Build the summary operation based on the TF collection of Summaries.
    summary_op = tf.merge_all_summaries()

    # Create a saver for writing training checkpoints.
    saver = tf.train.Saver()
    
    # Instantiate a SummaryWriter to output summaries and the Graph.
    summary_writer = tf.train.SummaryWriter(FLAGS.train_dir, sess.graph)

    # And then after everything is built, start the training loop.
    for step in xrange(FLAGS.max_steps):
      start_time = time.time()

      # Fill a feed dictionary with the actual set of images and targets
      # for this particular training step.
      feed_dict = fill_feed_dict(data_sets.train,
                                 feature_placeholder, 
                                 path_placeholder,
                                 targets_placeholder,
                                 sequence_length_placeholder,
                                 contextual_vector_placeholder)

      
      # Run one step of the model.  The return values are the activations
      # from the `train_op` (which is discarded) and the `loss` Op.  To
      # inspect the values of your Ops or variables, you may include them
      # in the list passed to sess.run() and the value tensors will be
      # returned in the tuple from the call.

      _, loss_value = sess.run([train_op, loss],
                               feed_dict=feed_dict)

      duration = time.time() - start_time

      # Write the summaries and print an overview fairly often.
      if step % 100 == 0:
        # Print status to stdout.
        print('Step %d: loss = %.2f (%.3f sec)' % (step, loss_value, duration))
        # Update the events file.
        summary_str = sess.run(summary_op, feed_dict=feed_dict)
        summary_writer.add_summary(summary_str, step)
        summary_writer.flush()

      # Save a checkpoint and evaluate the model periodically.
      if (step + 1) % 1000 == 0 or (step + 1) == FLAGS.max_steps:
        saver.save(sess, FLAGS.train_dir, global_step=step)
        # Evaluate against the training set.
        print('Training Data Eval:')
        do_eval(sess,
                eval_correct,
                feature_placeholder, 
                path_placeholder,
                targets_placeholder,
                sequence_length_placeholder,
                contextual_vector_placeholder,
                data_sets.train)
        # Evaluate against the validation set.
        print('Validation Data Eval:')
        do_eval(sess,
                eval_correct,
                feature_placeholder, 
                path_placeholder,
                targets_placeholder,
                sequence_length_placeholder,
                contextual_vector_placeholder,
                data_sets.validation)
        # Evaluate against the test set.
        print('Test Data Eval:')
        do_eval(sess,
                eval_correct,
                feature_placeholder, 
                path_placeholder,
                targets_placeholder,
                sequence_length_placeholder,
                contextual_vector_placeholder,
                data_sets.test)


def main(_):
  run_training()


if __name__ == '__main__':
  tf.app.run()
