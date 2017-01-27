from __future__ import print_function
import csv
import numpy as np
import argparse
from scipy.stats import pearsonr

def permute_data(data):
    data_len = len(data)
    perm = np.random.permutation(data_len)
    data_perm = data[perm]
    return data_perm


def cross_validation_split(data, crossval_split_index,
                           crossval_total_num_splits,
                           validation_data_ratio=0.1):
    assert validation_data_ratio > 0 and validation_data_ratio < 1
    assert crossval_split_index < crossval_total_num_splits

    N = len(data)
    n_test = int(N * 1. / crossval_total_num_splits)
    if crossval_split_index == crossval_total_num_splits - 1:
        n_test = N - crossval_split_index * n_test

    start_test = crossval_split_index * n_test
    end_test = crossval_split_index * n_test + n_test
    testdata = data[start_test: end_test]
    rest_data = np.concatenate((data[:start_test], data[end_test:]))

    n_valid = int(N * validation_data_ratio)
    valdata = rest_data[: n_valid]
    traindata = rest_data[n_valid:]
    return traindata, valdata, testdata

def read_csv(filename, smile_name, target_name, logp_name=None):
    data = []
    dt = np.dtype('S100, float')
    with open(filename) as file:
        reader = csv.DictReader(file)
        for row in reader:
            data_point=(row[smile_name], float(row[target_name]))
            if logp_name:
                dt = np.dtype('S100, float, float')
                data_point += (float(row[logp_name]),)
            data.append(data_point)

    return np.asarray(data, dtype=dt)


def model_params(s):
    try:
        x, y, z = map(int, s.split(','))
        return x, y, z
    except:
        raise argparse.ArgumentTypeError("Model paramaters must be x,y,z")


def get_metric(predictions, targets):
    rmse = np.sqrt(np.mean((predictions - targets) ** 2))
    aae = np.mean(np.abs(predictions - targets))
    r,_ = pearsonr(predictions,targets)
    return rmse, aae, r