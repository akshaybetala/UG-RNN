import numpy as np
import utils

# Load the data points D, and the station locations (lat/lon)
D = np.genfromtxt('data/data.txt',delimiter=None)
loc = np.genfromtxt('data/locations.txt',delimiter=None)
edges = np.genfromtxt('data/edges.txt',delimiter=None)
m,n = D.shape
# m = 2760 data points, n=30 dimensional
# D[i,j] = 1 if station j observed rainfall on day i

prob = np.mean(D,axis=0)
prob = [1-prob ,prob]

factors = utils.calculate_factors(D)
EMI_score = utils.calculate_emis(factors,prob)
max_spanning_tree = utils.maximum_spanning_tree(EMI_score*-1)
print "max_spanning_tree:" , max_spanning_tree
# utils.plot_graph(max_spanning_tree,loc)

log_likelihood = utils.caluculate_likelihood(D,factors,prob,max_spanning_tree)

print "log_likelihood:" , log_likelihood


