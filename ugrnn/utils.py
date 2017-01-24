from __future__ import print_function
import csv
import numpy as np


def permute_data(data, labels):
    data_len = len(data)
    perm = np.random.permutation(data_len)
    data_perm = data[perm]
    labels_perm = labels[perm]
    return data_perm, labels_perm


def cross_validation_split(data, labels, crossval_split_index,
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
    testdata = (data[start_test: end_test], labels[start_test: end_test])

    rest_data = np.concatenate((data[:start_test], data[end_test:]), axis=0)
    rest_labels = np.concatenate((labels[:start_test], labels[end_test:]), axis=0)

    n_valid = int(N * validation_data_ratio)
    valdata = (rest_data[: n_valid], rest_labels[: n_valid])
    traindata = (rest_data[n_valid:], rest_labels[n_valid:])
    print(len(traindata[0]), len(valdata[0]), len(testdata[0]))
    return traindata, valdata, testdata


def read_csv(filename, smile_name, target_name, logp_name=None):
    if logp_name:
        data = ([], [], [])
    else:
        data = ([], [])
    with open(filename) as file:
        reader = csv.DictReader(file)
        for row in reader:
            data[0].append(row[smile_name])
            data[1].append(float(row[target_name]))
            if logp_name:
                data[2].append(float(row[logp_name]))
    return map(np.array, data)




