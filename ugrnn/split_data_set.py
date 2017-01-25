import numpy as np

from utils import read_csv, cross_validation_split, permute_data

csv_file_path = 'ugrnn/data/delaney/delaney.csv'
smile_col_name = "smiles"
target_col_name = "solubility"
logp_col_name = "logp"

data = read_csv(csv_file_path, smile_col_name, target_col_name, logp_col_name)
data_perm = permute_data(data)

traindata, valdata, testdata = cross_validation_split(data, crossval_split_index=1, crossval_total_num_splits=10)

train_file_path = 'ugrnn/data/delaney/train_delaney.csv'
validate_file_path = 'ugrnn/data/delaney/validate_delaney.csv'
test_file_path = 'ugrnn/data/delaney/test_delaney.csv'


header = "{:},{:},{:}".format(smile_col_name, target_col_name, logp_col_name )
fmt = ('%s', '%4f', '%4f')
np.savetxt(train_file_path, traindata, header=header, fmt=fmt, comments='', delimiter=',')
np.savetxt(validate_file_path, valdata, header=header, fmt=fmt, comments='', delimiter=',')
np.savetxt(test_file_path, testdata, header=header, fmt=fmt, comments='', delimiter=',')

