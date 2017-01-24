from __future__ import print_function
import os
import numpy as np


def read_file(path):
	a =  np.loadtxt(path, usecols=(3,5,6), dtype={'names': ('label', 'logP', 'smiles'),
												 'formats': (np.float, np.float, 'S200')})
	print(a)


def read_data_set(add_logP=False):
	path = os.path.dirname(os.path.realpath(__file__))
	train_path = os.path.join(path, 'train.txt')
	test1_path = os.path.join(path, 'test1.smi')
	test2_path = os.path.join(path, 'test2.smi')

	train_data = read_file(train_path)
	test1_data = read_file(test1_path)
	test2_data = read_file(test2_path)


read_data_set()

