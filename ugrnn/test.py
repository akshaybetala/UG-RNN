# import tensorflow as tf
#
# with tf.Graph().as_default():
# 	sess = tf.Session()
# 	sequence_len = tf.constant(5)
#
#
# 	feature_pl  = tf.constant([[[1,1],[2,2],[3,3],[4,4], [5,5]],
# 								[[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]],
# 							    [[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]],
# 								[[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]],
# 								[[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]],])
#
# 	indices = tf.pack([tf.range(0, sequence_len), tf.range(0, sequence_len)], axis=1)
# 	step_feature = tf.squeeze(tf.gather_nd(feature_pl, indices))
#
# 	init = tf.global_variables_initializer()
# 	sess.run(init)
#
# 	x = sess.run([step_feature, indices])
# 	print x

import numpy as np



data_len = len(data)
perm = np.random.permutation(data_len)

