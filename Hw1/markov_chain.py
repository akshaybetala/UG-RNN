import numpy as np
from os import walk
mypath = 'proteins/'  # use path to data files
_, _, filenames = next(walk(mypath), (None, None, []))

np.set_printoptions(precision=2)

mSeq = len(filenames)        # read in each sequence
# mSeq = 10
o,x = [],[]
for i in range(mSeq):
    f = open('proteins/' +  filenames[i] , 'r')
    o.append( f.readline()[:-1] )  # strip trailing '\n'
    x.append( f.readline()[:-1] )
    f.close()

xvals, ovals = set(),set()  # extract the symbols used in x and o
for i in range(mSeq):
    xvals |= set(x[i])
    ovals |= set(o[i])
xvals = list( np.sort( list(xvals) ) )
ovals = list( np.sort( list(ovals) ) )
dx,do = len(xvals),len(ovals)

for i in range(mSeq):       # and convert to numeric indices
    x[i] = np.array([xvals.index(s) for s in x[i]])
    o[i] = np.array([ovals.index(s) for s in o[i]])

p0 = np.zeros(dx)
for i in range(mSeq):
    p0[x[i][0]] += 1

p0 = p0/sum(p0)
print 'p0:', p0

Tr = np.zeros((dx,dx))
for seq in range(mSeq):
    for s in range(len(x[seq])-1):
        Tr[x[seq][s]][x[seq][s+1]] += 1
Tr = Tr/Tr.sum(axis=1)[:,None]
print 'Tr:', Tr[:5,:5]

print np.matmul(p0,np.linalg.matrix_power(Tr,100))


Ob = np.zeros((dx,do))

for seq in range(mSeq):
    for s in range(len(x[seq])):
        Ob[x[seq][s]][o[seq][s]]+=1
Ob = np.asarray(Ob/Ob.sum(axis=1)[:,None])

print 'Ob:', Ob[:5,:5]

# o = [1,2,3]

def markovMarginals(o,p0,Tr,Ob):
    '''Compute p(o) and the marginal probabilities p(x_t|o) for a Markov model
       defined by P[xt=j|xt-1=i] = Tr(i,j) and P[ot=k|xt=i] = Ob(i,k) as numpy matrices'''
    dx,do = Ob.shape   # if a numpy matrix
    L = len(o)
    f = np.zeros((L,dx))
    r = np.zeros((L,dx))
    p = np.zeros((L,dx))
    f[0,:] = p0*Ob[:,o[0]]    # compute initial forward message
    log_pO =  np.log(f[0,:].sum())  # update probability of sequence so far
    f[0,:] /= f[0,:].sum()  # normalize (to match definition of f)

    for t in range(1,L):    # 	compute forward messages
    	f[t,:] = np.matmul(f[t-1,:],Tr)*Ob[:,o[t]]
        log_pO += np.log(f[t,:].sum())
        f[t,:] /= f[t,:].sum()
        
    r[L-1,:] = np.ones(dx)  # initialize reverse messages
    p[L-1,:] = r[L-1,:]*f[L-1,:]  # and marginals

  
    for t in range(L-2,-1,-1):
        r[t,:] =  np.matmul(Tr,r[t+1,:]*Ob[:,o[t+1]])
        r[t,:] /= r[t,:].sum()
        p[t,:] = r[t,:]*f[t,:]
        p[t,:] /= p[t,:].sum()

    return log_pO, p

def testMarkovMarginals():
	Tr = np.asarray([[0, 0, 1],[.33, .66, 0], [.5, .5, 0]])
	Ob = np.asarray([[1,0],[.5,.5],[0,1]])
	p0 = np.asarray([.33, .33, .33])
	O = [0,1]
	log_pO, p = markovMarginals(O,p0, Tr,Ob)
	print p
	
testMarkovMarginals()
_, p = markovMarginals(o[0],p0,Tr,Ob)
print p[6,:]

_, p = markovMarginals(o[2],p0,Tr,Ob)
print p[9,:]

log_pO, _ = markovMarginals([1,2,3,4,5],p0,Tr,Ob)
print log_pO


print filenames[0]
print filenames[2]
print filenames[4]