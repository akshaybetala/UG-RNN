import tensorflow as tf

x = tf.ones([1000,100,10], dtype=tf.float32, name=None)
x.set_shape([None,None,10])
y = tf.split(0,1, x)
z = tf.zeros([100,10],dtype=tf.float32)

for i,j in enumerate(y):
	print i
	z = tf.add(j,z)

# Create a session for running Ops on the Graph.
sess = tf.Session()
y_value = z.eval(session=sess)
print y_value.shape()