sudo apt-get update
sudo apt-get install python-pip python-dev
sudo apt-get install python-rdkit librdkit1 rdkit-data
sudo pip install --upgrade networkx 
export TF_BINARY_URL=https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-0.11.0rc0-cp27-none-linux_x86_64.whl
sudo pip install --upgrade $TF_BINARY_URL
