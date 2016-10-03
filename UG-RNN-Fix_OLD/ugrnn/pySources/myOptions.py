class myOptions:
	def __init__(self,optFilePath,models):
		self.models = models
		self.classes=0
		self.ONN_OUTPUTS_NUMBER=0
		self.ONN_HIDDEN_UNITS=[]
		try:
			self.optFile = open(optFilePath,"r")
		except IOError:
			exit("CANNOT OPEN myOptions FILE")
		self.read()
		self.optFile.close()
	def read(self):
		optFileList = self.optFile.readlines()
		self.classes = optFileList[0].split(" ")[1]
		if int(self.classes) != 1:
			exit("classes MUST BE EQUAL TO 1")
		self.ONN_OUTPUTS_NUMBER = optFileList[1].split(" ")[1]
		try:
			if not(int(self.ONN_OUTPUTS_NUMBER) > 0):
				self.optFile.close()
				exit("ONN_OUTPUTS_NUMBER MUST BE GREATER THAN 0")
		except Exception:
			self.optFile.close()
			exit("\nCHECK ONN_OUTPUTS_NUMBER IN myOptions FILE")	 
		self.ONN_HIDDEN_UNITS = optFileList[2].split(" ")[1:self.models+1]
		try:
			check = [int(x) > 0 for x in self.ONN_HIDDEN_UNITS]
			for truth in check:
				if not truth:
					self.optFile.close()
					exit("ONN_HIDDEN_UNITS MUST BE GREATER THAN 0")
		except Exception:
			self.optFile.close()
			exit("\nCHECK ONN_HIDDEN_UNITS IN myOptions FILE")
		self.MNN_OUTPUTS_NUMBER= optFileList[3].split(" ")[1:self.models+1]
		self.MNN_HIDDEN_UNITS = optFileList[4].split(" ")[1:self.models+1] 
		self.seed = optFileList[5].split(" ")[1:self.models+1]
		self.epsilon = optFileList[6].split(" ")[1:self.models+1]
		self.nEpochs = optFileList[7].split(" ")[1]
		self.batchBlocks= optFileList[8].split(" ")[1]
		self.folds= optFileList[9].split(" ")[1]
		self.epsilonMin = optFileList[10].split(" ")[1:self.models+1]
		self.epsEpoch = optFileList[11].split(" ")[1]
		self.properties = optFileList[12].split(" ")[1]
	def write(self,filePath,model):
		try:
			outFile = open(filePath,"w")
		except:
			exit("CANNOT OPEN FILE "+outFile)
		outFile.write("classes "+self.classes)
		outFile.write("ONN_OUTPUTS_NUMBER "+self.ONN_OUTPUTS_NUMBER)
		if((model+1)<self.models):
			outFile.write("ONN_HIDDEN_UNITS "+self.ONN_HIDDEN_UNITS[model]+"\n")
		else:	
			outFile.write("ONN_HIDDEN_UNITS "+(self.ONN_HIDDEN_UNITS[model]).strip("\n")+"\n")
		if((model+1)<self.models):
			outFile.write("MNN_OUTPUTS_NUMBER "+self.MNN_OUTPUTS_NUMBER[model]+"\n")
		else:	
			outFile.write("MNN_OUTPUTS_NUMBER "+(self.MNN_OUTPUTS_NUMBER[model]).strip("\n")+"\n")
		if((model+1)<self.models):
			outFile.write("MNN_HIDDEN_UNITS "+self.MNN_HIDDEN_UNITS[model]+"\n")
		else:	
			outFile.write("MNN_HIDDEN_UNITS "+(self.MNN_HIDDEN_UNITS[model]).strip("\n")+"\n")
		if((model+1)<self.models):
			outFile.write("seed "+self.seed[model]+"\n")
		else:	
			outFile.write("seed "+(self.seed[model]).strip("\n")+"\n")
		if((model+1)<self.models):
			outFile.write("epsilon "+self.epsilon[model]+"\n")
		else:	
			outFile.write("epsilon "+(self.epsilon[model]).strip("\n")+"\n")
		outFile.write("nEpochs "+self.nEpochs+"\n")
		outFile.write("batchBlocks "+self.batchBlocks)
		outFile.write("folds "+self.folds)
		if((model+1)<self.models):
			outFile.write("epsilonMin "+self.epsilonMin[model]+"\n")
		else:	
			outFile.write("epsilonMin "+(self.epsilonMin[model]).strip("\n")+"\n")
		outFile.write("epsEpoch "+self.epsEpoch)
		outFile.write("properties "+self.properties)
		outFile.close()
