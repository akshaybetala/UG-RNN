import numpy as np
import pyGM as gm
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix

def calculate_factors(D):
	m,n = D.shape
	# define binary random variables
	X = [gm.Var(i,2) for i in range(n)]

	# create empty factor over binary x0,x1
	factors = [[gm.Factor([X[i],X[j]],0.0) for j in xrange(0,n)] for i in xrange(0,n)]
	for i in xrange(0,n-1):
		for j in xrange(i+1,n):
			if i==j:
				continue
			for s in D:
				factors[i][j].table[s[i],s[j]]  +=1

	for i in xrange(0,n):
		for j in xrange(0,n):
			if i!=j:
				factors[i][j].table /= factors[i][j].table.sum()

	return factors

def calculate_emis(factors,prob):
	n = len(prob[0])
	EMI_score = np.zeros((n,n))

	def emis_of_edge(x1, x2, factors, prob):
		score = 0
		for i in [0,1]:
			for j in [0,1]:
				join_probability = factors[x1][x2].table[i,j]
				score+= join_probability*np.log(join_probability/(prob[i][x1]*prob[j][x2]))
		return score

	for i in xrange(0,n-1):
		for j in xrange(i+1,n):
			EMI_score[i,j] = EMI_score[j,i] =  emis_of_edge(i,j,factors,prob)
	return EMI_score

def caluculate_likelihood(D,factors,prob,edges):
	log_likelihood = 0
	m,n = D.shape

	def joint_likelihood(x1, x2, v1, v2, factors):
		x1 = int(x1)
		x2 = int(x2)
		v1 = int(v1)
		v2 = int(v2)
		join_probability = factors[x1][x2].table[v1,v2]
		return np.log(join_probability/(prob[v1][x1]*prob[v2][x2]))


	for s in D:
		log_likelihood = 0
		for i in xrange(0,n):
			log_likelihood +=np.log(prob[int(s[i])][i])
		for edge in edges:
			x1 = int(edge[0])
			x2 = int(edge[1])
			log_likelihood+= joint_likelihood(x1,x2,s[x1],s[x2],factors)
	return log_likelihood/m

def plot_graph(edges,loc):

	x_cord = loc[:,0]
	y_cord = loc[:,1]

	for edge  in edges:
		plt.plot( [x_cord[edge[0]], x_cord[edge[1]]] , [y_cord[edge[0]], y_cord[edge[1]]],marker='o' )

	plt.show()

def maximum_spanning_tree(weight):
	n = len(weight) 
	X = csr_matrix(weight)
	Tcsr = minimum_spanning_tree(X)
	temp = Tcsr.toarray()
	final_edges = []
	for i in xrange(n):
		for j in xrange(n):
			if(temp[i,j]<0):
				if i<j:
					final_edges.append([i,j])
				else:
					final_edges.append([j,i])
	return final_edges
