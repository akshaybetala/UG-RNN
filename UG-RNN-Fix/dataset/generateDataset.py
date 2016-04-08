#Auhtor: Alessandro Lusci
import openbabel
import pybel
import tools 
import sys
import os

def getRingAtomsMulti(molfilePath):
	for m in pybel.readfile("mol", molfilePath):
		#print m.OBMol.GetFormula()
		outString = ""
		rings = m.OBMol.GetSSSR()
		numR = 0;
		for r in rings:
			numR = numR+1
		outString+="Ring count "+str(numR)+"\n"
		for r in rings:
			outString+="Ring size "+str(r.Size())+"\n"
			path = r._path
			for p in path:
				outString+=str(p-1)+"\n"

		outString+="Atom count "+str(m.OBMol.NumAtoms())+"\n"
		outString+="Index HowManyRings RingSize Hybridization Hydro_count Aromaticity AntiClockwise_chiral"+"\n"
		i = 0
		for a in openbabel.OBMolAtomIter(m.OBMol):
			outString+=str(i)+" "+str(a.MemberOfRingCount())+" "+str(a.MemberOfRingSize())+" "+str(a.GetHyb())+" "+str(a.ImplicitHydrogenCount())+" "
			if a.IsAromatic():
				outString+="1"+" "
			else:
				outString+="0"+" "
			if a.IsClockwise():
				outString+="1"+"\n"
			else:
				outString+="0"+"\n"
			i = i + 1

	return outString


#main
usageString = "USAGE: python generateDatset.py -d dataset -r"
if (len(sys.argv)<3 or sys.argv[1]!="-d"):
	exit(usageString)
datasetName = sys.argv[2]
reduceGraph = 0
addLogP = 1
for i in range(3,len(sys.argv)):
	if sys.argv[i] == "-r":
		reduceGraph = 1
#open files
try:
	smileFile = open(datasetName+"/"+datasetName+".smi","r")
except IOError:
	exit("CANNOT OPEN SMILE FILE: "+datasetName+"/"+datasetName+".smi")
try:
	targetFile = open(datasetName+"/"+datasetName+".target","r")
except IOError:
	exit("CANNOT OPEN TARGET FILE: "+datasetName+"/"+datasetName+".target")
try:
	logPFile = open(datasetName+"/"+datasetName+".logP","r")
except IOError:
	addLogP = 0
	print("NO LOGP FILE FOUND: "+datasetName+"/"+datasetName+".logP")
try:
	datasetFile = open(datasetName+"/"+datasetName+"Temp","w")
except IOError:
	exit("CANNOT OPEN DATASET FILE: "+datasetName+"/"+datasetName+"Temp")
###########
#count molecules in dataset
smileList = smileFile.readlines()
targetList = targetFile.readlines()
if addLogP == 1:
	logPList = logPFile.readlines()
totalMolecules = len(smileList)
if totalMolecules != len(targetList):
	exit("ERROR: smileFile and targetFile contain a different number of molecules")
if addLogP == 1 and totalMolecules!=len(logPList):
	exit("ERROR: smileFile and logPFile contain a different number of molecules")
datasetFile.write(str(totalMolecules)+"\n")
###########################
#main loop
for l in range(0,len(smileList)):
	print "molecule "+str(l)
	print smileList[l]
	smile = smileList[l].strip("\n")
	target = targetList[l].strip("\n")
	logP = "0.0"	
	if addLogP == 1:
		logP = logPList[l].strip("\n")
	molecule = pybel.readstring("can",smile)
	#count atoms number in molecule
	atomsNumber = 0
	for atom in molecule:
		atomsNumber+=1
	if atomsNumber >= 100:
		print "Warning molecule "+str(l)+" contains "+str(atomsNumber)+" atoms"
	###############################
	datasetFile.write("train"+str(l)+".can")
	molecule.write("mol",datasetName+"/temp","overwrite=true")	
	
	if reduceGraph==0:
		tools.standardizeMol(datasetName)		
		try:
			temp = open(datasetName+"/temp","r")
		except IOError:
			exit("CANNOT OPEN TEMP FILE")
		#delete M RAD and M CHG lines
		tools.deleteLines(datasetName)
		#standardize molfile
		tools.standardizeMol(datasetName)		
		datasetFile.write(temp.read()+"$$$$\n"+target+"\n")		
		datasetFile.write(logP+"\n")
				
	else:
		#try:
		#	temp = open("temp","r")
		#except IOError:
		#	exit("CANNOT OPEN TEMP FILE")
		#standardize molfile		
		ringsInfo = getRingAtomsMulti(datasetName+"/temp")
		print ringsInfo
		#standardize molfile
		tools.standardizeMol(datasetName)		
		#temp.close()
		try:
			temp = open(datasetName+"/temp","a")
		except IOError:
			exit("CANNOT OPEN TEMP FILE")
		temp.write("$$$$\n"+target+"\n"+ringsInfo)
		temp.close()
		
		#delete M RAD and M CHG lines
		tools.deleteLines(datasetName)
		#reduce graph here
		tools.reduceGraph(datasetName)
		
		try:
			tempReduced = open(datasetName+"/tempReduced","r")
		except IOError:
			exit("CANNOT OPEN "+datasetName+"/tempReduced FILE")
		datasetFile.write(tempReduced.read())
		datasetFile.write(logP+"\n")
		tempReduced.close()
		os.remove(datasetName+"/tempReduced")		

os.remove(datasetName+"/temp")
datasetFile.close()
##########
#####
