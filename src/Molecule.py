import numpy as np 
from rdkit import Chem
from rdkit.Chem import Draw
import parse_solubility_data
import networkx as nx

class Molecule:
	graph = None
	atoms = []
	no_of_atoms = 0
	directed_graphs = {}
	m = None

	def __init__(this,smile):
		this.m = Chem.MolFromSmiles(smile)
		Chem.Kekulize(this.m)
		this.no_of_atoms = this.m.GetNumAtoms()
		this.graph = nx.Graph()
		for i in xrange(this.no_of_atoms):
			this.graph.add_node(i)
			atom = this.m.GetAtomWithIdx(i)
			this.atoms.append(atom.GetSymbol())
			for neighbour in atom.GetNeighbors():
				neighbour_idx = neighbour.GetIdx()
				this.graph.add_edge(i,neighbour_idx)

	def get_bond_type_as_double(this,idx1,idx2):
		bond = this.m.GetBondBetweenAtoms(i,neighbour_idx)
		if bond==None:
			return 0
		else:
			return bond.GetBondTypeAsDouble()

	def create_directed_graphs(this):
		for idx in xrange(this.no_of_atoms):
			path = nx.single_source_dijkstra_path(this.graph,idx)
			G = nx.DiGraph()
			for i in xrange(this.no_of_atoms):
				temp = path[i]
				temp.reverse()
				G.add_path(temp)
			this.directed_graphs[idx] = nx.topological_sort(G)
