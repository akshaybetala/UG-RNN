#Author Alessandro Lusci
import sys
import os
import math
import random
from pySources.train import *
from pySources.myOptions import myOptions as myOpt 


#main
def main(argv):
	#parse input
	usageString = "python ugrnn.py -d dataset -t testID -m numberOfUGRNNs -f numberOfFolds \noptional parameters:\n-submit \n-train [0:train and validate,1:train validate and test,2:train validate and test with fixed sequences,default=0] \n-vp % (default=20%)"
	if (len(argv)<9):
		exit(usageString)
	
	#parameters
	train = Train(argv,usageString) 
	
	#get number of molecules in dataset
	train.getMoleculesInDataset()	
	
	#generate directories structure
	train.generateDirectoriesStructure()
	
	#write option files
	train.writeOptFiles()
	
	#generate molecules sequences
	train.generateSequences()
	
	#write sequences
	train.writeSequences()

	#create PBS dir
	train.createPBSDir()

	#submit jobs
	if train.submit == 1:
		train.submitJobsToQueue()
		
		
if __name__=="__main__":
	main(sys.argv)
