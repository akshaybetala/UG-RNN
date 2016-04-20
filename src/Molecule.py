import numpy as np 
from rdkit.Chem.rdchem import BondType
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
		self.create_initial_vectors()

	def create_directed_graphs(self):
		self.directed_graphs = {}
		# self.bonds_list = {}

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
			sorted_path = []
			# bonds = []
			for i in xrange(self.no_of_atoms-1):
				node = topological_order[i]
				sorted_path.append((nx.edges(G,node))[0])
				# bonds.append(get_bond_type_as_double())

			self.directed_graphs[idx] = sorted_path[:]

	def create_initial_vectors(self):
		# create a three dimesnional matrix G, such that Gij is the contextual vector for ith vertex in jth DAG
		self.initial_vector = np.zeros((self.no_of_atoms,self.no_of_atoms,71))
		for idx in xrange(self.no_of_atoms):
			sorted_path = self.directed_graphs[idx]
			for i,j in sorted_path:
				# print np.append(self.get_atom_type_as_one_hot_vector(i) ,self.get_bond_type_as_one_hot_vector(i,j))
				self.initial_vector[idx,i,:] = np.append(self.get_atom_type_as_one_hot_vector(i) ,self.get_bond_type_as_one_hot_vector(i,j))

	def get_bond_type_as_one_hot_vector(self,idx1,idx2):
		bond_type = self.get_bond_type(idx1,idx2)
		one_hot = np.zeros(21)
		one_hot[bond_type]=1
		return one_hot

	def get_bond_type(self,idx1,idx2):
		bond = self.m.GetBondBetweenAtoms(idx1,idx2)
		bond_type = bond.GetBondType()

		if(bond_type == BondType.UNSPECIFIED):
			return 0
		elif(bond_type == BondType.SINGLE):
			return 1
		elif(bond_type == BondType.DOUBLE):
			return 2
		elif(bond_type == BondType.TRIPLE):
			return 3
		elif(bond_type == BondType.QUADRUPLE):
			return 4
		elif(bond_type == BondType.QUINTUPLE):
			return 5
		elif(bond_type == BondType.HEXTUPLE):
			return 6
		elif(bond_type == BondType.ONEANDAHALF):
			return 7
		elif(bond_type == BondType.TWOANDAHALF):
			return 8
		elif(bond_type == BondType.THREEANDAHALF):
			return 9
		elif(bond_type == BondType.FOURANDAHALF):
			return 10
		elif(bond_type == BondType.FIVEANDAHALF):
			return 11
		elif(bond_type == BondType.AROMATIC):
			return 12
		elif(bond_type == BondType.IONIC):
			return 13
		elif(bond_type == BondType.HYDROGEN):
			return 14
		elif(bond_type == BondType.THREECENTER):
			return 15
		elif(bond_type == BondType.DATIVEONE):
			return 16
		elif(bond_type == BondType.DATIVE):
			return 17
		elif(bond_type == BondType.DATIVEL):
			return 18
		elif(bond_type == BondType.DATIVER):
			return 19
		elif(bond_type == BondType.OTHER):
			return 20

	def get_atom_type_as_one_hot_vector(self,idx):
		one_hot = np.zeros(50)
		atom_type = self.get_atom_type(idx)
		one_hot[atom_type]=1
		return one_hot

	def get_atom_type(self,idx):
		s = self.atoms[idx]

		if (s == 'C'):
			return 0
		elif (s == 'H'):
			return 1
		elif (s == 'N'):
			return 2 
		elif (s == 'O'):
			return 3 
		elif (s == 'S'):
			return 4 
		elif (s == 'F'):
			return 5 
		elif (s == 'Cl'):
			return 6 
		elif (s == 'Br'):
			return 7
		elif (s == 'I'):
			return 8 
		elif (s == 'In'):
			return 9 			
		elif (s == 'K'):
			return 10
		elif (s == 'Na'):
			return 11
		elif (s == 'Ba'):
			return 12
		elif (s == 'Sb'):
			return 13
		elif (s == 'P'):
			return 14	
		elif (s == 'Be'):
			return 15	
		elif (s == 'Sn'):
			return 16
		elif (s == 'Cu'):
			return 17
		elif (s == 'B'):
			return 18
		elif (s == 'Cd'):
			return 19
		elif (s == 'Ca'):
			return 20
		elif (s == 'As'):
			return 21
		elif (s == 'Co'):
			return 22
		elif (s == 'Cr'):
			return 23
		elif (s == 'Te'):
			return 24
		elif (s == 'Fe'):
			return 25
		elif (s == 'Pb'):
			return 26
		elif (s == 'Mn'):
			return 27
		elif (s == 'Hg'):
			return 28
		elif (s == 'Mo'):
			return 29
		elif (s == 'Ni'):
			return 30
		elif (s == 'Se'):
			return 31
		elif (s == 'Ti'):
			return 32
		elif (s == 'Zn'):
			return 33
		elif (s == 'Si'):
			return 34
		elif (s == '?'):
			return 35
		elif (s == 'Mg'):
			return 36
		elif (s == 'V'):
			return 37 
		elif (s == 'Li'):
			return 38 
		elif (s == 'Al'):
			return 39
		elif (s == 'Zr'):
			return 40
		elif (s == 'Bi'):
			return 41
		elif (s == 'Pd'):
			return 42
		elif (s == 'Pt'):
			return 43
		elif (s == 'Ru'):
			return 44
		elif (s == 'Rh'):
			return 45			
		elif (s == 'Ga'):
			return 46						
		elif (s == 'Ge'):
			return 47								
		elif (s == 'Ag'):
			return 48
		# elif (s == 'Tb'):
			# return 49
		# elif (s == 'Ir'):
		# 	return 50
		# elif (s == 'W'):
		# 	return 51
		# elif (s == 'Cs'):
		# 	return 52
		# elif (s == 'Re'):
		# 	return 53
		# elif (s == 'Pr'):
		# 	return 54
		# elif (s == 'Nd'):
		# 	return 55
		# elif (s == 'Gd'):
		# 	return 56
		# elif (s == 'Yb'):
		# 	return 57
		# elif (s == 'Er'):
		# 	return 58
		# elif (s == 'U'):
		# 	return 59
		# elif (s == 'Tl'):
		# 	return 60
		# elif (s == 'Au'):
		# 	return 61
		# elif (s == 'Ac'):
		# 	return 62
		# elif (s == 'Ho'):
		# 	return 63
		# elif (s == 'Os'):
		# 	return 64
		# elif (s == 'Sm'):
		# 	return 65
		# elif (s == 'Nb'):
		# 	return 66
		else:
			return 49		

	