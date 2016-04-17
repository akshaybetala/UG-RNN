import numpy as np 
from rdkit import Chem
from rdkit.Chem import Draw
import parse_solubility_data
import networkx as nx

class Molecule:
	
	def __init__(self,smile,size):
		self.atoms = []
		self.m = Chem.MolFromSmiles(smile)
		Chem.Kekulize(self.m)
		self.no_of_atoms = self.m.GetNumAtoms()
		self.graph = nx.Graph()
		for i in xrange(self.no_of_atoms):
			self.graph.add_node(i)
			atom = self.m.GetAtomWithIdx(i)
			self.atoms.append(atom.GetSymbol())
			for neighbour in atom.GetNeighbors():
				neighbour_idx = neighbour.GetIdx()
				self.graph.add_edge(i,neighbour_idx)

		self.create_directed_graphs()
		self.create_contextual_vectors(size)


	def get_bond_type_as_double(self,idx1,idx2):
		bond = self.m.GetBondBetweenAtoms(i,neighbour_idx)
		if bond==None:
			return 0
		else:
			return bond.GetBondTypeAsDouble()

	def create_directed_graphs(self):
		self.directed_graphs = {}
		for idx in xrange(self.no_of_atoms):
			path = nx.single_source_dijkstra_path(self.graph,idx)
			G = nx.DiGraph()
			for i in xrange(self.no_of_atoms):
				temp = path[i]
				temp.reverse()
				G.add_path(temp)
			self.directed_graphs[idx] = nx.topological_sort(G)

	def create_contextual_vectors(self,vector_length):
		# create a three dimesnional matrix G, such that Gij is the contextual vector for ith vertex in jth DAG
		self.contextual_vector = np.zeros((self.no_of_atoms,self.no_of_atoms,vector_length))