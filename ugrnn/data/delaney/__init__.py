from __future__ import print_function
import os
import numpy as np

def generate_random_permute():
	path = os.path.dirname(os.path.realpath(__file__))
	smiles_path = os.path.join(path, 'delaney.smi')
	target_path = os.path.join(path, 'delaney.target')
	logp_path = os.path.join(path, 'delaney.logP')

	with open(smiles_path) as f:
		smiles = np.array([line.rstrip() for line in f])

	with open(target_path) as f:
		labels = np.array([float(line.rstrip()) for line in f])

	with open(logp_path) as f:
		logP = np.array([float(line.rstrip()) for line in f])

	data_len = len(smiles)
	perm = np.random.permutation(data_len)

	perm_smiles = smiles[perm]
	perm_labels = labels[perm]
	perm_logP = logP[perm]

	perm_smiles_path = os.path.join(path, 'delaney_perm.smi')
	perm_target_path = os.path.join(path, 'delaney_perm.target')
	perm_logp_path = os.path.join(path, 'delaney_perm.logP')

	np.savetxt(perm_smiles_path,perm_smiles, fmt="%s")
	np.savetxt(perm_target_path, perm_labels, fmt="%f")
	np.savetxt(perm_logp_path, perm_logP, fmt="%f")

def read_data_set(add_logP=False):
	path = os.path.dirname(os.path.realpath(__file__))
	smiles_path = os.path.join(path, 'delaney.smi')
	target_path = os.path.join(path, 'delaney.target')
	logp_path = os.path.join(path, 'delaney.logP')

	perm_smiles_path = os.path.join(path, 'delaney_perm.smi')
	perm_target_path = os.path.join(path, 'delaney_perm.target')
	perm_logp_path = os.path.join(path, 'delaney_perm.logP')

	with open(perm_smiles_path) as f:
		smiles = np.array([line.rstrip() for line in f])

	with open(perm_target_path) as f:
		labels = np.array([float(line.rstrip()) for line in f])

	logP = None

	if add_logP:
		with open(perm_logp_path) as f:
			logP = np.array([float(line.rstrip()) for line in f])

	return smiles, labels, logP

