import sys
import os
from myOptions import * 
from utility import *
import math
import random

class Train:
	datasetName = ""
	testID = ""
	models = 0
	folds = 0
	train = 0
	submit = 0
	valPerc = 20
	testSeq = []
	trainSeq = []
	validationSeq = []
	
	def __init__(self,argv,usageString):
		for i in range(0,len(argv)):
			if argv[i] == "-d":
				self.datasetName = sys.argv[i+1]
			elif argv[i] == "-t":
				self.testID = sys.argv[i+1]
			elif argv[i] == "-m":
				self.models = int(sys.argv[i+1])
			elif argv[i] == "-f":	
				self.folds = int(sys.argv[i+1])
			elif argv[i] == "-submit":
				self.submit = 1			
			elif argv[i] == "-train":
				self.train = int(argv[i+1])			
			elif argv[i] == "-vp" and validation == 1:
				self.valPerc = int(argv[i+1])
				print "VALIDATION_PERC == "+argv[i+1]			
				if valPerc < 5 or valPerc > 35:
					exit("VALIDATION_PERC OUT OF BOUND: IT MUST BE BETWEEN 5 AND 35")						  			  			
		if (self.datasetName == "" or self.testID == "" or self.models <= 0 or self.folds <= 0):
			exit(usageString)
						
	def getMoleculesInDataset(self):
		try:
			self.datasetFile = open("./data/"+self.datasetName,"r")
		except IOError:
			exit("CANNOT OPEN DATASET "+self.datasetName)
		#get number of molecules
		self.molecules = int(self.datasetFile.readline())
		print "UGRNN TRAIN ON "+self.datasetName+" DATASET\nNUMBER OF MOLECULES "+str(self.molecules)
		self.datasetFile.close()
	
	def generateDirectoriesStructure(self):		
		print "GENERATING DIRECTORY STRUCTURE..."
		if os.system("mkdir test/"+self.testID):
			exit("DIRECTORY test/"+self.testID+" ALREADY EXISTS\nRUN rm -r test/"+self.testID+" TO GO ON")
		for f in range(0,self.folds):
			os.system("mkdir test/"+self.testID+"/fold"+str(f))
			for m in range(0,self.models):
				os.system("mkdir test/"+self.testID+"/fold"+str(f)+"/model"+str(m))
		print "...DONE"
	
	def writeOptFiles(self):
		print "WRITING OPT FILES..."
		self.myOptions = myOptions("myOptions",self.models)
		if int(self.myOptions.folds.strip("\n")) != self.folds:
			exit("NUMBER OF FOLD IN OPTION FILE DOES NOT MATCH THE NUMBER YOU ENTERED!"+"\nmyOptions.folds="+str(self.myOptions.folds.strip("\n"))+" folds="+str(self.folds))
		for m in range(0,self.models):
			self.myOptions.write("test/"+self.testID+"/fold0"+"/model"+str(m)+"/myOptions",m)
		print "...DONE"			

	def trainAndValidate(self):		
		for f in range(0,self.folds):
			self.validationSeq.append([])
		m = 0
		while m < self.molecules:
			for f in range(0,self.folds):
				if m < self.molecules:
					self.validationSeq[f].append(str(m))
				m+=1			
		for f in range(0,self.folds):
			self.trainSeq.append([])
		m = 0
		while m < self.molecules:
			for f in range(0,self.folds):
				if find(m,self.validationSeq[f]) == 0:
					self.trainSeq[f].append(str(m))
			m += 1
		
	def trainValidateAndTest(self):
		trainSeqTemp = []
		for f in range(0,self.folds):
			self.testSeq.append([])
		m = 0
		while m < self.molecules:
			for f in range(0,self.folds):
				if m < self.molecules:
					self.testSeq[f].append(str(m))
				m+=1
		for f in range(0,self.folds):
			trainSeqTemp.append([])
		m = 0
		while m < self.molecules:
			for f in range(0,self.folds):
				if find(m,self.testSeq[f]) == 0:
					trainSeqTemp[f].append(str(m))
			m += 1
		for f in range(0,self.folds):
			self.trainSeq.append([])
			self.validationSeq.append([])
		for f in range(0,self.folds):			
			step = int(len(trainSeqTemp[f])/int(self.valPerc*len(trainSeqTemp[f])/100))			
			for m in range(0,len(trainSeqTemp[f])):
				if math.fmod(m,step)==0:
					self.validationSeq[f].append(str(trainSeqTemp[f][m]))
				else:
					self.trainSeq[f].append(str(trainSeqTemp[f][m]))
	
	def trainValidatAndTestFixedSeq(self):
		trainLength = 0
		validationLength = 0
		testLength = 0
		seed = 3
		try:
			sequencesFile = open("sequences")
		except IOError:	
			exit("CANNOT OPEN sequences file")
		sequences = sequencesFile.readlines()
 		trainLength = int((sequences[0].split(" ")[1]).strip("\n"))
		testLength = int((sequences[1].split(" ")[1]).strip("\n"))
		trainSeqTemp = []
		for f in range(0,self.folds):
			self.trainSeq.append([])
			trainSeqTemp.append([])
			self.validationSeq.append([])
			self.testSeq.append([])
			for m in range(0,trainLength):
				trainSeqTemp[f].append(str(m))
			random.seed(f+seed)
			random.shuffle(trainSeqTemp[f])		
			step = int(len(trainSeqTemp[f])/int(self.valPerc*len(trainSeqTemp[f])/100))
			for m in range(0,len(trainSeqTemp[f])):
				if math.fmod(m,step)==0:
					self.validationSeq[f].append(str(trainSeqTemp[f][m]))
				else:
					self.trainSeq[f].append(str(trainSeqTemp[f][m]))
			for m in range(trainLength,trainLength+testLength):
				self.testSeq[f].append(str(m))
			
	def generateSequences(self):
		if self.train == 0:
			self.trainAndValidate()
		elif self.train == 1:
			self.trainValidateAndTest()
		elif self.train == 2:
			self.trainValidatAndTestFixedSeq()	
			
	def writeSequences(self):
		for f in range(0,self.folds):
			if self.train == 0:
				print "fold"+str(f)+"\nTRAINING SET LENGTH:"+str(len(self.trainSeq[f]))+"\nVALIDATION SET LENGTH:"+str(len(self.validationSeq[f]))
			else:
				print "fold"+str(f)+"\nTRAINING SET LENGTH:"+str(len(self.trainSeq[f]))+"\nVALIDATION SET LENGTH:"+str(len(self.validationSeq[f]))+"\nTEST SET LENGTH:"+str(len(self.testSeq[f]))
		
		for f in range(0,self.folds):
			try:
				trainSeqFile = open("test/"+self.testID+"/fold"+str(f)+"/train","w")
			except IOError:
				exit("CANNOT OPEN test/"+self.testID+"/fold"+str(f)+"/train")
			trainSeqFile.write(str(len(self.trainSeq[f]))+"\n");
			trainSeqFile.write(" ".join(self.trainSeq[f]))
			trainSeqFile.close()
			try:
				validationFile = open("test/"+self.testID+"/fold"+str(f)+"/validation","w")
			except IOError:
				exit("CANNOT OPEN test/"+self.testID+"/fold"+str(f)+"/validation")
			validationFile.write(str(len(self.validationSeq[f]))+"\n");
			validationFile.write(" ".join(self.validationSeq[f]))
			validationFile.close()
			if self.train != 0:			
				try:
					testSeqFile = open("test/"+self.testID+"/fold"+str(f)+"/test","w")
				except IOError:
					exit("CANNOT OPEN test/"+self.testID+"/fold"+str(f)+"/test")
				testSeqFile.write(str(len(self.testSeq[f]))+"\n");
				testSeqFile.write(" ".join(self.testSeq[f]))
				testSeqFile.close()
	
	def createPBSDir(self):
		os.system("mkdir test/"+self.testID+"/PBS")
		for f in range(0,self.folds):
			for m in range(0,self.models):
				try:
					qsubScript = open("test/"+self.testID+"/PBS/"+self.testID+str(f)+"F"+str(m)+"M","w")
				except IOError:
					exit("CANNOT OPEN PBS FILE")
				qsubScript.write("#!/bin/bash\n")
				qsubScript.write("#PBS -l walltime=100:00:00\n")
				qsubScript.write("#PBS -e error.txt\n")
				qsubScript.write("#PBS -o output.txt\n")
				qsubScript.write("#PBS -V\n")
				qsubScript.write("cd "+os.getcwd()+"\n")
				if self.train == 0:
					qsubScript.write("./ugrnnTrain -d "+self.datasetName+" -t "+self.testID+" -f "+str(f)+" -m "+str(m)+" -vl > ./test/"+self.testID+"/fold"+str(f)+"/model"+str(m)+"/logModel"+str(m))									
				else:
					qsubScript.write("./ugrnnTrain -d "+self.datasetName+" -t "+self.testID+" -f "+str(f)+" -m "+str(m)+" -tsval > ./test/"+self.testID+"/fold"+str(f)+"/model"+str(m)+"/logModel"+str(m))				
				qsubScript.close()
				
	def submitJobsToQueue(self):
		#run tests
		for f in range(0,self.folds):
			for m in range(0,self.models):
				os.system("qsub ./test/"+self.testID+"/PBS/"+self.testID+str(f)+"F"+str(m)+"M")	
