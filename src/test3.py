import tensorflow as tf

def myfunction(vector):

  result = tf.reduce_sum(vector)
  # print_result = tf.Print(result, [result], "myfunction called ")
  return result

MAX_ROWS = 10

# input matrix with 2 columns and unknown number of rows (<MAX_ROWS)
inputs = tf.placeholder(tf.int32, [None, 2])
# copy of inputs, will need to have a persistent copy of it because we will
# be fetching rows in different session.run calls
# data = tf.Variable(inputs, validate_shape=False)
# input producer that iterates over the rows and pushes them onto Queue
rows = tf.train.slice_input_producer([inputs], shuffle=False)
result = tf.zeros([0],tf.int32)
for row in rows:
	result = tf.add(tf.reduce_sum(row),result)

myfunction_op = result

# this op will save placeholder values into the variable
init_op = tf.initialize_all_variables()

# Coordinator is not necessary in this case, but you'll need it if you have
# more than one Queue in order to close all queues together
sess = tf.Session()
coord = tf.train.Coordinator()
threads = tf.train.start_queue_runners(sess=sess, coord=coord)

sess.run([init_op],feed_dict={inputs:[[0, 0], [1, 1], [2, 2]]})

try:
  # for i in range(MAX_ROWS):
  a = sess.run([myfunction_op])
  print(a)
except tf.errors.OutOfRangeError:
  print('Done iterating')
finally:
  # When done, ask other threads to stop.
  coord.request_stop()