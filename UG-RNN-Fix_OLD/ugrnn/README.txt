AUTHOR: ALESSANDRO LUSCI UGRNN-R SOFTWARE PACKAGE
VERSION 1.0
Software: UG-RNN SOFTWARE PACKAGE Copyright(c) 2013 All rights reserved

The software included in this release are provided as is without implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See LICENSE
file included in this release.

Ubuntu is a registered trademark of Canonical Ltd.  
python-openbabel is part of the Open Babel project
CentOS is part of the Community ENTerprise Operating System
ChemAxon Products include tools for visualization and drawing of molecules, 
chemical database searching and management, drug discovery. 
ChemAxon provides products free for academic use.

--SOFTWARE REQUIREMENTS

Ubuntu 9.10 or higher 
CentOS 5.8 or higher
python-openbabel (http://openbabel.org/wiki/Python)

-----------------------

--COMPILE THE SOURCE CODE

cd ugrnn/ 
make 
make -f makefilePredict 

--HOW TO TRAIN A MODEL

In order to train a model you need a dataset of molecules in the UGRNN format.
The UGRNN format is a molfile V2000 format, with the addition of the desidered
target and the logP. In order to convert a dataset of molecules to the UGRNN
format, you can run the following bash script:

/bin/bash generateDataset.sh datasetName -r

datasetName is a directory that must contain two files: datasetName.smi and
datasetName.target. The first one is a list of molecular smiles.The second
one is the list of related targets. The script generates a dataset which is
named datasetName and it is placed under the directory datasetName. If logP values
are available and you want to train a model with logPs, put a file contaning the
list of LogP values under datasetName directory. The file must be named
datasetName.LogP and logP values has to be in the same order as in
datasetName.smi.
-r option generates molfiles where the rings are contracted to single nodes.

In order to compute logPs you can use many softwares from the literature. 
We chose ChemAxon Marvin beans from http://www.chemaxon.com/products/marvin/ 
  
After the dataset has been generated, you have to move it into the directory
/ugrnn/data. 

Now you are ready to run the script:

python ugrnnTrain.py -d dataset -t testID -m numberOfUGRNNs -f numberOfFolds
optional parameters: -submit -train [0:train and validate,1:train validate and
test,2:train validate and test with fixed sequences,default=0] -vp %
(default=20%)     
      
ugrnnTrain.py requires four parameters: a dataset file that must be under
ugrnn/data directory; a testID which is the name of the directory where the
training results are saved; the number of ugrnn models that you want to be part
of the NN-ensemble; number of folds for cross-validation (if this number is
equal to the number of molecules in the traning set, a leave-one-out
cross-validation is performed). We recommend to start with a 10-fold cross
validation. By default, the script perform a training procedure on 80% of the
total dataset and a validation procedure on the remaining 20%. Both training set
and validation set are randomly selected from the original dataset.  If you want
to assess the performance of you trained model after on an external test set,
select the option train 1. The program will randomly split the original dataset
into a training set and a test and 20% of the traning set will be randomly
selected for validation purposes.  
If you wish to assess the perfomances of the trained models, use the
option -train 1.  
If you wish to test the same sequence for each fold, use
option -train 2. 
This option requires a sequence file, that has to be placed
under sequence directory and has to be in the following format:

train 60 
test 14
        
Morover, under ugrnn directory, it is necessary to create a myOptions file which
specifies the Neural Networks parameters.  
myOptions must be in the following format:

classes 1 
ONN_OUTPUTS_NUMBER 1 
ONN_HIDDEN_UNITS 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 
MNN_OUTPUTS_NUMBER 3 4 5 6 7 8 9 10 11 12 3 3 3 3 3 3 3 3 3 3
MNN_HIDDEN_UNITS 7 7 7 7 7 7 7 7 7 7 3 4 5 6 7 8 9 10 11 12 
seed 10 11 12 13 14 15 16 17 18 19 20 10 10 10 10 10 10 10 10 10 
epsilon 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01
nEpochs 5000 
batchBlocks 25 
folds 10 
epsilonMin 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 
epsEpoch 5000

Notice that each number in myOptions file is separated by a space. 
batchBlocks specifies the number of molecule per batch training. 

In the current version of the software classes can only be 1. 
Classes is a flag that has been added to allow classification predictions but 
the current UG-RNN version has been tested only on regression tasks. 
      
----------------------

--HOW TO CHECK RESULTS

python ugrnnCheckResults.py

Usage: python ugrnnCheckResults.py -t testID -m numberOfModels -f numberOfFolds
-vl -aof

testID is the absolute path to the directory which contains the trained models.
m is the number of models that have been trained.  
f is the number of folds into which the dataset has been split.  
select option -aof when you want to compute the average over folds
 
Printed results refer to the following metrics: 
Absoulte Average Error (AAE)
Root Mean Square Error (RMSE) 
Squared Pearson Correlation Coefficient (PEARSON)
 
----------------------

--HOW TO MAKE PREDICTIONS

python ugrnnPredict.py

Usage: ugrnnPredict [OPTIONS]
		-i	input molecules
		-md	model directory
		-m	number of models						
		-f	number of folds
		-r	enable graph reduction

In order to run the script ugrnnPredict.py you need: 
a file containing a list of smiles (type the file name after "-i"); 
one or more trained models (type the model directory after "-md"); 
In addition to this  you have to specify the number of trained models ("-m" option) 
per fold and the total number of folds ("-f" option) 
into which the training set has been split.  
"-r" is an option that enables graph reduction and must be used 
only if models have been trained on a dataset where graphs have been reduced. 
Don not use "-r" otherwise.
   
-------------------------

--EXAMPLE

Generate a UGR-NN-dataset 
######################## 
cd UGRNN-R/dataset/ 
/bin/bash generateDataset.sh delaney 
mv delaney/delaney ../ugrnn/data

Train a model on the delaney dataset 
##################################### 
cd UGRNN-R/ugrnn 
python ugrnnTrain.py -d delaney -t delaney -m 1 -f 10

---OUTPUT---

UGRNN TRAIN AND TEST ON delaney DATASET
NUMBER OF MOLECULES 1144
GENERATING DIRECTORY STRUCTURE...
...DONE
WRITING OPT FILES...
...DONE
fold0
TRAINING SET LENGTH:1029
VALIDATION SET LENGTH:115
fold1
TRAINING SET LENGTH:1029
VALIDATION SET LENGTH:115
fold2
TRAINING SET LENGTH:1029
VALIDATION SET LENGTH:115
fold3
TRAINING SET LENGTH:1029
VALIDATION SET LENGTH:115
fold4
TRAINING SET LENGTH:1030
VALIDATION SET LENGTH:114
fold5
TRAINING SET LENGTH:1030
VALIDATION SET LENGTH:114
fold6
TRAINING SET LENGTH:1030
VALIDATION SET LENGTH:114
fold7
TRAINING SET LENGTH:1030
VALIDATION SET LENGTH:114
fold8
TRAINING SET LENGTH:1030
VALIDATION SET LENGTH:114
fold9
TRAINING SET LENGTH:1030
VALIDATION SET LENGTH:114

---OUTPUT---

The above script creates all files and directories which are necessary to train
a model but it does not start the actual training, unless the -submit option is
selected. If so, the script tries to submit 10 jobs to a torque queuing system.
Obviously, the script will fail if torque is not installed or not properly configured. 
For information on how to install torque, check out the following page:
http://wiki.debian.org/Torque    

Anyway you can always train a model on a specific fold, running the command
string that is contanined in the pbs file
UG-RNN/ugrnn/test/delaney/PBS/testDelaneyfFmM.  

For example, testDelaney0F0M will look like this:
#!/bin/bash 
#PBS -l walltime=100:00:00 
#PBS -e error.txt 
#PBS -o output.txt 
#PBS -V 
cd /home/yourHomeDir/UGRNN-R/ugrnn 
./ugrnnTrain -d delaney -t delaney -f 0 -m 0 -vl > ./test/delaney/fold0/model0/logModel0 

If you run the command:
nohup ./ugrnnTrain -d delaney -t delaney -f 0 -m 0 -vl > ./test/delaney/fold0/model0/logModel0 &

a training process of model 0 on fold 0 starts.  

Check the training logs using the following command:
tail -f ./test/delaney/fold0/model0/logModel0

In order train the model on each fold, you have to iterate over each fold.

./ugrnnTrain -d delaney -t testDelaney -f 1 -m 0 -vl 
./ugrnnTrain -d delaney -t testDelaney -f 2 -m 0 -vl 
./ugrnnTrain -d delaney -t testDelaney -f 3 -m 0 -vl
./ugrnnTrain -d delaney -t testDelaney -f 4 -m 0 -vl 
./ugrnnTrain -d delaney -t testDelaney -f 5 -m 0 -vl 
./ugrnnTrain -d delaney -t testDelaney -f 6 -m 0 -vl
./ugrnnTrain -d delaney -t testDelaney -f 8 -m 0 -vl 
./ugrnnTrain -d delaney -t testDelaney -f 9 -m 0 -vl 
./ugrnnTrain -d delaney -t testDelaney -f 10 -m 0 -vl

Notice that the above commands can run in parallel.
##################################### 

Check results on validation set
###############################
In oreder to check the results on the validation set during the training process,
run the following command:
python ugrnnCheckResults.py -t test/delaney -f 1 -m 1 -vl -aof
Completed folds: 1
Fold0
RMSE:1.12305501671 AAE:0.921557638261 PEARSON:0.676731964261
FOLD AVERAGE...
RMSE:1.12305501671 0.0 AAE:0.921557638261 0.0 PEARSON:0.676731964261 0.0

Make a prediction on a trained model 
#################################### 
When the training process ends, you are able to predict logS values of external
molecules (i.e. molecules that are not included in the training set).  
For example, if you want to predict the log solubility values of the molecules
contained in the file input.smi, you can run the following command:

python ugrnnPredict.py -i input.smi -md test/delaney -m 1 -f 1

--OUTPUT--

molecule 0
ClCC(Cl)(Cl)Cl
molecule 1
CC(Cl)(Cl)Cl
molecule 2
ClC(Cl)C(Cl)Cl
molecule 3
ClCC(Cl)Cl
molecule 4
FC(F)(Cl)C(F)(Cl)Cl
molecule 5
CC(Cl)Cl
molecule 6
ClC(=C)Cl
molecule 7
CCOC(C)OCC
molecule 8
Clc1ccc(Cl)c(Cl)c1Cl
molecule 9
C1CCc2ccccc2C1
PREDICTING...
RESULTS
fold0
molecule 0 -2.83102
molecule 1 -2.1448
molecule 2 -3.22657
molecule 3 -2.51806
molecule 4 -2.37556
molecule 5 -1.83564
molecule 6 -1.7054
molecule 7 -1.41659
molecule 8 -4.79288
molecule 9 -3.72698
AVERAGE RESULTS
molecule 0 -2.83102
molecule 1 -2.1448
molecule 2 -3.22657
molecule 3 -2.51806
molecule 4 -2.37556
molecule 5 -1.83564
molecule 6 -1.7054
molecule 7 -1.41659
molecule 8 -4.79288
molecule 9 -3.72698

--OUTPUT--

The above program computes the solubility value of each molecule, using a model
trained on fold 0. Of course, if you trained a model on every fold, you could run
the command:
python ugrnnPredict.py -i inputSmile -md test/testDelaney -m 1 -f 10  

The output of the above command will be the average, over 10 folds, of the
solubility values of each molecule.  
####################################

---------
