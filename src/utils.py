import numpy as np 
from rdkit.Chem.rdchem import BondType
from rdkit import Chem

def num_of_features():
		return num_bond_features()+num_atom_features()

def atom_features(symbol):
	return np.array(one_of_k_encoding_unk(symbol,
		['C','H','N', 'O', 'S', 'F','Cl','Br','I','In','K','Na',
		'Ba','Sb','P','Be','Sn','Cu','B','Cd','Ca','As','Co','Cr',
		'Te','Fe','Pb','Mn','Hg','Mo','Ni','Se','Ti','Zn','Si',
		'Mg','V','Li','Al','Zr','Bi','Pd','Pt','Ru','Rh','Ga','Ge',
		'Ag','Tb','Ir','W','Cs','Re','Pr','Nd','Gd','Yb','Er',
		'U','Tl','Au','Ac','Ho','Os','Sm','Nb','?']))

				# one_of_k_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5]) +
				# one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4]) +
				# one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5]) +
				# [atom.GetIsAromatic()])

def bond_features(bond):
	bt = bond.GetBondType()
	return np.array([bt == BondType.UNSPECIFIED,
					bt == BondType.SINGLE,
					bt == BondType.DOUBLE,
					bt == BondType.TRIPLE,
					bt == BondType.QUADRUPLE,
					bt == BondType.QUINTUPLE,
					bt == BondType.HEXTUPLE,
					bt == BondType.ONEANDAHALF,
					bt == BondType.TWOANDAHALF,
					bt == BondType.THREEANDAHALF,
					bt == BondType.FOURANDAHALF,
					bt == BondType.FIVEANDAHALF,
					bt == BondType.AROMATIC,
					bt == BondType.IONIC,
					bt == BondType.HYDROGEN,
					bt == BondType.THREECENTER,
					bt == BondType.DATIVEONE,
					bt == BondType.DATIVE,
					bt == BondType.DATIVEL,
					bt == BondType.DATIVER,
					bt == BondType.OTHER])

def one_of_k_encoding(x, allowable_set):
	if x not in allowable_set:
		raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
	return map(lambda s: x == s, allowable_set)

def one_of_k_encoding_unk(x, allowable_set):
	"""Maps inputs not in the allowable set to the last element."""
	if x not in allowable_set:
		x = allowable_set[-1]
	return map(lambda s: x == s, allowable_set)

def num_atom_features():
	# Return length of feature vector using a very simple molecule.
	m = Chem.MolFromSmiles('CC')
	alist = m.GetAtoms()
	a = alist[0]
	return len(atom_features(a))

def num_bond_features():
	# Return length of feature vector using a very simple molecule.
	simple_mol = Chem.MolFromSmiles('CC')
	Chem.SanitizeMol(simple_mol)
	return len(bond_features(simple_mol.GetBonds()[0]))