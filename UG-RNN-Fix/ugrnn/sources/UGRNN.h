#ifndef UGRNN_H 
#define UGRNN_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "Molecule.h"
#include "NNt.h"
#include "NNr.h"
#include "myOptions.h"
#include "General.h"
#include "Node.h"
#include "nodeIO.h"
#include <cmath>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <time.h>

class UGRNN{
	protected:
	NNt* ONN;
	NNt* MNN;
	int ATOM_NUMBER;
	int MAX_CHILDREN;
	int TOTAL_INPUTS;
	int propertiesNum;
	int firstEpoch;
	vector<string> totalAtoms;
	Float averageAtomsNumber;
	myOptions Opt;
	vector<int> trainVector;
	vector<int> testVector;	
	vector<int> validationVector;
	vector<int> inputVector;
	string modelDir;
	string logFilePath; 
	ofstream logTP;
	ofstream logTPVal;
	ofstream logError;
	ofstream logTestError;
	ofstream logValidationError;
	ofstream logTestVector;
	ofstream logTrainDescriptors;
	ofstream logTestDescriptors;
	public:
	//empty constructor
	UGRNN();
	//constructor
	UGRNN(NNr* _ONN,NNt* _MNN,Molecule** moleculeList);
	//constructor_2
   	UGRNN(Molecule** moleculeList,myOptions  & Opt,int moleculeNumber,string modelDir ,int predictor=0);
   	//copy constructor
	UGRNN(const UGRNN & _UGRNN);
	//destructor
	~UGRNN();	
	//assignment operator
	UGRNN operator=(const UGRNN & _uGRNN);
	//getters
	NNt* getONN(){return ONN;}
	NNt* getMNN(){return MNN;}
	int getATOM_NUMBER(){return ATOM_NUMBER;}
	int getMAX_CHILDREN(){return MAX_CHILDREN;}
	int getTOTAL_INPUTS(){return TOTAL_INPUTS;}
	int getPropertiesNum(){return propertiesNum;}
	int getFirstEpoch(){return firstEpoch;}
	vector<string> getTotalAtoms(){return totalAtoms;}
	Float getAverageAtomsNumber(){return averageAtomsNumber;}
	vector<int> getTrainVector(){return trainVector;}
	vector<int> getTestVector(){return testVector;}
	vector<int> getValidationVector(){return validationVector;}
	vector<int> getInputVector(){return inputVector;}
	string getLogFilePath(){return logFilePath;}
	string getModelDir(){return modelDir;}
	ofstream & getLogTP(){return logTP;}
	ofstream & getLogTPVal(){return logTPVal;}
	ofstream & getLogError(){return logError;}
	ofstream & getLogTestError(){return logTestError;}
	ofstream & getLogValidationError(){return logValidationError;}
	ofstream & getLogTestVector(){return logTestVector;}
	ofstream & getLogTrainDescriptors(){return logTrainDescriptors;}
	ofstream & getLogTestDescriptors(){return logTestDescriptors;}
	//setters
	void setONN(NNt* _ONN){ONN = _ONN;}
	void setMNN(NNt* _MNN){MNN = _MNN;}
	void setATOM_NUMBER(int _ATOM_NUMBER){ATOM_NUMBER=_ATOM_NUMBER;}
	void setMAX_CHILDREN(int _MAX_CHILDREN){MAX_CHILDREN=_MAX_CHILDREN;}
	void setTOTAL_INPUTS(int _TOTAL_INPUTS){TOTAL_INPUTS=_TOTAL_INPUTS;}
	void setPropertiesNum(int _propertiesNum){propertiesNum=_propertiesNum;}
	void setFirstEpoch(int _firstEpoch){firstEpoch=_firstEpoch;}
	void setTotalAtoms(vector<string> _totalAtoms){totalAtoms=_totalAtoms;}
	void setTrainVector(vector<int> _trainVector){trainVector = _trainVector;}
	void setTestVector(vector<int> _testVector){testVector = _testVector;}
	void setValidationVector(vector<int> _validationVector){validationVector = _validationVector;}
	void setInputVector(vector<int> _inputVector){inputVector = _inputVector;}
	void setLogFilePath(string _logFilePath){logFilePath = _logFilePath;}
	void setModelDir(string _modelDir){modelDir = _modelDir;}
	void setLogTP(string logTPPath){
		logTP.open(logTPPath.c_str(),ios_base::app);
		//debug
		cout << logTPPath << endl;
		///////
		if(!logTP.is_open()){
			cerr << "CANNOT OPEN logTP file" << endl;
			exit(1);
		}
	}
	void setLogTPVal(string logTPValPath){
		logTPVal.open(logTPValPath.c_str(),ios_base::app);
		//debug
		cout << logTPValPath << endl;
		///////
		if(!logTPVal.is_open()){
			cerr << "CANNOT OPEN logTPVal file" << endl;
			exit(1);
		}
	}
	void setLogError(string logErrorPath){
		logError.open(logErrorPath.c_str(),ios_base::app);
		if(!logError.is_open()){
			cerr << "CANNOT OPEN logError file" << endl;
			exit(1);
		}
	}
	void setLogTestError(string logTestErrorPath){
		logTestError.open(logTestErrorPath.c_str(),ios_base::app);
		if(!logTestError.is_open()){
			cerr << "CANNOT OPEN logTestError file" << endl;
			exit(1);
		}
	}
	void setLogValidationError(string logValidationErrorPath){
		logValidationError.open(logValidationErrorPath.c_str(),ios_base::app);
		if(!logValidationError.is_open()){
			cerr << "CANNOT OPEN logTestError file" << endl;
			exit(1);
		}
	}
	void setLogTestVector(string logTestVectorPath){
		logTestVector.open(logTestVectorPath.c_str(),ios_base::app);
		if(!logTestVector.is_open()){
			cerr << "CANNOT OPEN logTestVector file" << endl;
			exit(1);
		}
	}
	void setLogTrainDescriptors(string logTrainDescriptorsPath){
		logTrainDescriptors.open(logTrainDescriptorsPath.c_str());
		if(!logTrainDescriptors.is_open()){
			cerr << "CANNOT OPEN trainDescriptors file" << endl;
			exit(1);
		}
	}
	void setLogTestDescriptors(string logTestDescriptorsPath){
		logTestDescriptors.open(logTestDescriptorsPath.c_str());
		if(!logTestDescriptors.is_open()){
			cerr << "CANNOT OPEN testDescriptors file" << endl;
			exit(1);
		}
	}
	//untranslateU: translate IAN's atom symbol code into its original string name
	string unTranslateU(int a);
	//getAllAtoms: return a vector containing all atoms in the dataset, identified by their chemical name
	vector<string> getAllAtoms(Molecule** moleculeList,int moleculeNumber);  
	//transform an atom symbol into its bit representation
	int* atomToArray(string atomSymbol);
	//create model dir
	void createModelDir(string  & modelDir,int fold);
	//read model: takes ONNMoldePath and MNNModelPath
	void read(string ONNModelPath,string MNNModelPath);
	//write model: takes ONN and MNN ostream
	void write(ofstream  & model,ofstream & modelM);	
	//train and validate
	void trainAndValidate(Molecule** train,Molecule** validation,const vector<Float> & bound,int numberOfModels);
	//validate
	void validate(Molecule** validate,const vector<Float> & bound,Float* lowestAEE,vector<Float> & bestP,vector<Float> & bestT,int e, int numberOfModelsi,Float** AAE);
	//predict
	void predict(Molecule** testSet, const vector<Float> & bound,int numberOfModels);
	//predict target
	vector<Float> predictTarget(Molecule** molecules, const vector<Float> bound, int numberOfModels);
	//forward through MNN and ONN
	Float forward(Molecule* molecule,Float* properties,Float* descriptors=NULL);
	//compute: computes the output for a node in the graph
	Float* compute(int node,NodeIO** nodeIO,int** children,int* num_children,int numChildren,int** bond_type);
	//compute back: performs gradient descent
	void computeBack(int node,NodeIO** nodeIO,Float* delta,int** children,int* num_children,int numChildren); 
	//initNodeIO
	void initNodeIO(Molecule* molecule,NodeIO*** nodeIO);
	//forwardMNN
	void forwardMNN(Molecule* molecule,NodeIO*** nodeIO,Float* out);
	//backwardMNN
	void backwardMNN(Molecule* molecule,NodeIO*** nodeIO);	
	//initialize categorical input
	int getMaxNumberOfChildren(Molecule** moleculeList,int moleculeNumber);
	//accumulate outSubGraph
	void sum(Float* out,Float* outSubgraph,int outputNumber);
	//compute the average number of atoms per molecule
	Float computeAverageAtomsNumber(Molecule** moleculeArray,int moleculeNumber); 
	//normalize MNN OUTPUTS
	void normalizeOut(Float* out,int outputNumber,int totalAtoms);
	//compute average out
	void computeAverageOut(Float* totalOut,int outputNumber,int moleculeNumber);
	//shuffle inputs
	static vector<int> shuffle(int moleculeNumber,int seed);
	//init weights
	void initWeights(int seed);	
	//print float array
	void print(Float* array,int outputNumber,string message);
	//print train and test descriptors to files
	//FILES: trainDescriptor1... 
	void printDescriptors(Molecule** train,Molecule** test,vector<Float> bound,int numberOfModels);	
	//compute correlation coefficient r = (sigma(1,N)ti*pi-N*average(t)*average(p))/((N-1)*st*sp)
	//st = sqrt(sigma(N,1)(ti-average(t))^2)
	//sp = sqrt(sigma(N,1)(pi-average(p))^2)
	//this is a Pearson's correlation coefficient
	static Float computeCorrelation(const vector<Float> & t,const vector<Float> & p,int N);
	//compute RMSE and AAE
	static vector<Float> computeError(const vector<Float> & t,const vector<Float> & p,int N);		
	//denormalizing function
	static Float denormalize(Float value,const vector<Float> & bound);
	//normalize logP
	Float normalizeLogP(Float value,const vector<Float> & bound);
	//write totalAtoms
	void writeTotalAtoms(string dir);
	//read totalAtoms
	void readTotalAtoms(string dir);
	//compare totalAtoms
	void compareTotalAtoms(vector<string> totalAtoms,vector<string> totalAtomsTest);
	//write MAX_CHILDREN
	void writeMAX_CHILDREN(string dir);
	//read MAX_CHILDREN
	void readMAX_CHILDREN(string dir);
	
};
#endif
