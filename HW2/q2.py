import numpy as np
import utils

# Load the data points D, and the station locations (lat/lon)
D = np.genfromtxt('data/data.txt',delimiter=None)
loc = np.genfromtxt('data/locations.txt',delimiter=None)
edges = np.genfromtxt('data/edges.txt',delimiter=None)
m,n = D.shape
# m = 2760 data points, n=30 dimensional
# D[i,j] = 1 if station j observed rainfall on day i

# utils.plot_graph(edges)

prob = np.mean(D,axis=0)
prob = [1-prob ,prob]

factors = utils.calculate_factors(D)
EMI_score = utils.calculate_emis(factors,prob)

emi_score_of_edges = np.zeros(len(edges))
for i,edge in enumerate(edges):
	emi_score_of_edges[i] = EMI_score[edge[0],edge[1]]

print emi_score_of_edges