from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np


def load_solubility_data(file_path):
    f = file(file_path,'r')
    txt = f.read()
    f.close()
    #print txt
    lines = [line.split('"') for line in txt.split('\n')]
    parsed = []
    for line in lines[:-1]:
        if len(line)>1:
            parsed.append([x.replace(' ','') for x in ([line[1]]+line[2].split(',')) if x != ''])
        else:
            parsed.append([x.replace(' ','') for x in line[0].split(',')])
        
    prediction_targets = np.array([float(x[1]) for x in parsed[1:]],dtype='float32')
    SMILES = [x[-1] for x in parsed[1:]]
    return SMILES, prediction_targets
