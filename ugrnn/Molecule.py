from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np 
from rdkit.Chem.rdchem import BondType
from rdkit import Chem
from rdkit.Chem import Draw
import networkx as nx
import utils


class Molecule:

	num_of_features = 0

	def __init__(self,smile):
		self.smile = smile
		self.num_of_features = utils.num_of_features()
		self.atoms = []
		self.m = Chem.MolFromSmiles(smile)
		Chem.Kekulize(self.m)
		self.no_of_atoms = self.m.GetNumAtoms()
		self.graph = nx.Graph()
		
		for i in xrange(self.no_of_atoms):
			self.graph.add_node(i)
			atom = self.m.GetAtomWithIdx(i)
			print(atom.GetSymbol())
			self.atoms.append(atom)
			for neighbour in atom.GetNeighbors():
				neighbour_idx = neighbour.GetIdx()
				self.graph.add_edge(i,neighbour_idx)
		
		self.create_directed_graphs()
		self.create_feature_vectors()

	def create_directed_graphs(self):
		self.directed_graphs = np.empty((self.no_of_atoms,self.no_of_atoms-1,3),dtype=int)
	
		#parse all the atoms one by one and get directed graph to that atom
		for idx in xrange(self.no_of_atoms):
			#get shortest path from the root to all the other atoms and then reverse the edges.
			path = nx.single_source_dijkstra_path(self.graph,idx)
			G = nx.DiGraph()
			for i in xrange(self.no_of_atoms):
				temp = path[i]
				temp.reverse()
				G.add_path(temp)

			# do a topological sort to get a order of atoms with all edges pointing to the root
			topological_order = nx.topological_sort(G)
			sorted_path = np.empty((self.no_of_atoms-1,3))
			no_of_incoming_edges = {}
			for i in xrange(self.no_of_atoms-1):
				node = topological_order[i]
				edge = (nx.edges(G,node))[0]
				if edge[1] in no_of_incoming_edges:
					index = no_of_incoming_edges[edge[1]]
					no_of_incoming_edges[edge[1]] += 1
				else:
					index = 0
					no_of_incoming_edges[edge[1]] = 1
				sorted_path[i,:] = [node,edge[1],index]
			# sorted_path[self.no_of_atoms-1, :] = [idx,self.no_of_atoms]
			self.directed_graphs[idx,:,:] = sorted_path

	def create_feature_vectors(self):
		# create a three dimesnional matrix G, such that Gij is the contextual vector for ith vertex in jth DAG
		self.feature_vector = np.zeros((self.no_of_atoms,self.no_of_atoms,self.num_of_features))
		
		for idx in xrange(self.no_of_atoms):
			sorted_path = self.directed_graphs[idx,:,:]
			for i in xrange(self.no_of_atoms-1):
				node1 = sorted_path[i,0]
				node2 = sorted_path[i,1]
				self.feature_vector[idx,i,:] = np.append(utils.atom_features(self.atoms[node1]),
						utils.bond_features(self.m.GetBondBetweenAtoms(node1,node2)))

			self.feature_vector[idx, self.no_of_atoms-1,:] = np.append(utils.atom_features(self.atoms[idx]),
						np.zeros(utils.num_bond_features()))
	
