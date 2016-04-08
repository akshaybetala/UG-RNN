import os
import sys
import subprocess
from pySources.dataset import * 


def parseInput(argv):
	args = {"inputFile":"","modelDir":"","models":0,"folds":0,"reduceGraph":0}
	for i in range(0,len(argv)):
		if argv[i] == "-i":
			try:
				args["inputFile"] = argv[i+1]
			except IndexError,e:
				exit(e)
		elif argv[i] == "-md":
			try:
				args["modelDir"] = argv[i+1]
			except IndexError,e:
				exit(e)
		elif argv[i] == "-m":
			try:
				args["models"] = int(argv[i+1])
			except Exception,e:
				exit(e)
		elif argv[i] == "-f":
			try:				
				args["folds"] = int(argv[i+1])
			except Exception,e:
				exit(e)
		elif argv[i] == "-r":
			args["reduceGraph"] = 1
	return args
	
def predict(inputFile,modelDir,models,folds,molecules):	
	results = []
	for f in range(0,folds):	
		results.append([])	
		for r in range(0,molecules):
			results[f].append(0.0)
	for f in range(0,folds):
		for m in range(0,models):
			command = "".join(["./ugrnnPredict"," -i ",inputFile," -md ",modelDir," -f ",str(f)," -m ",str(m)])
			#print command
			proc = subprocess.Popen(["./ugrnnPredict","-i",inputFile,"-md",modelDir,"-f",str(f),"-m",str(m)],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdout,stderr = proc.communicate()
			if stderr:
				exit(stderr)
			else:
				out = stdout.split("\n")
				temp = out[len(out)-3].split(" ")
				for r in range(0,molecules):
					results[f][r] += float(temp[r])				
		for r in range(0,molecules):
			results[f][r] = results[f][r]/float(models)
	
	return results					
			
def main(argv):
	#parse input	
	"""Usage: ugrnnPredict [OPTIONS]
		-i	input molecules	[it must a file containing a list of smiles]
		-md	model directory
		-m	number of models						
		-f	number of folds
		-r	enable graph reduction"""
	usageString = """Usage: ugrnnPredict [OPTIONS]
		-i	input molecules
		-md	model directory
		-m	number of models						
		-f	number of folds
		-r	enable graph reduction"""
	args = parseInput(argv)
	if len(args.keys()) == 0 or args["models"] < 1 or args["folds"] < 1:
		exit(usageString)
	
	#generate dataset 
	molecules = generateInput(args)
	
	#predict property
	print "PREDICTING..."
	results = predict(args["inputFile"]+".mol",args["modelDir"],args["models"],args["folds"],molecules)
		
	#compute average results
	print("RESULTS")	
	avResults = []
	for r in range(0,molecules):
		avResults.append(0.0)	
	for f in range(0,args["folds"]):
		print "fold"+str(f)
		for r in range(0,molecules):
			print "molecule "+str(r)+" "+str(results[f][r])
			avResults[r] += results[f][r]
	avResults = [x/float(args["folds"]) for x in avResults]
	print "AVERAGE RESULTS"	
	for r in range(0,molecules):
		print "molecule "+str(r)+" "+str(avResults[r]) 
	os.system("rm temp")
		 	

if __name__ == "__main__":
	main(sys.argv)
