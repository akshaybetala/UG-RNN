#ifndef MYOPTIONS_H
#define MYOPTIONS_H
#include<iostream>
#include<fstream>
#include<sstream> 
#include<cstring>
#include<string> 
#include"General.h"
#include<stdlib.h>
using namespace std;
class myOptions{
	protected:
	int classes;
	int ONN_OUTPUTS_NUMBER;
	int ONN_HIDDEN_UNITS;
	int MNN_OUTPUTS_NUMBER;
	int MNN_HIDDEN_UNITS;
	int seed;
	Float epsilon;
	int nEpochs;
	int batchBlocks;
	int folds;
	Float epsilonMin;
	int epsEpoch;
    int properties;
	public:
	//empty constructor
	myOptions(){};
	//constructor to be updated
	myOptions(int _classes,int _ONN_OUTPUTS_NUMBER,int _ONN_HIDDEN_UNITS,int _MNN_OUTPUTS_NUMBER,int _MNN_HIDDEN_UNITS);
	//constructor
    myOptions(ifstream& is);
    //copy constructor
	myOptions(const myOptions & _Opt);
	//destructor
	~myOptions();
	//assignment operator
	myOptions& operator=(const myOptions& _Opt){
		classes = _Opt.classes;
        	ONN_OUTPUTS_NUMBER = _Opt.ONN_OUTPUTS_NUMBER;
       		ONN_HIDDEN_UNITS = _Opt.ONN_HIDDEN_UNITS;
        	MNN_OUTPUTS_NUMBER = _Opt.MNN_OUTPUTS_NUMBER;
        	MNN_HIDDEN_UNITS = _Opt.MNN_HIDDEN_UNITS;
        	seed =  _Opt.seed;
        	epsilon =  _Opt.epsilon;
        	nEpochs = _Opt.nEpochs;
        	batchBlocks =  _Opt.batchBlocks;
        	folds =  _Opt.folds;
        	epsilonMin = _Opt.epsilonMin;
			epsEpoch = _Opt.epsEpoch;
			properties = _Opt.properties;
	}
	//getters
	int getClasses(){return classes;}
	int getONN_OUTPUTS_NUMBER(){return ONN_OUTPUTS_NUMBER;}
	int getONN_HIDDEN_UNITS(){return ONN_HIDDEN_UNITS;}
	int getMNN_OUTPUTS_NUMBER(){return MNN_OUTPUTS_NUMBER;}
	int getMNN_HIDDEN_UNITS(){return MNN_HIDDEN_UNITS;}
	int getSeed(){return seed;}
	Float getEpsilon(){return epsilon;}
	int getNEpochs(){return nEpochs;}
	int getBatchBlocks(){return batchBlocks;}
	int getFolds(){return folds;}
	Float getEpsilonMin(){return epsilonMin;}
	int getEpsEpoch(){return epsEpoch;}
	int getProperties(){return properties;}
	//setters
	//to be updated
	void setClasses(int _classes){classes = _classes;}
	void setONN_OUTPUTS_NUMBER(int _ONN_OUTPUTS_NUMBER){ONN_OUTPUTS_NUMBER = _ONN_OUTPUTS_NUMBER;}
	void setONN_HIDDEN_UNITS(int _ONN_HIDDEN_UNITS){ONN_HIDDEN_UNITS = _ONN_HIDDEN_UNITS;}
	void setMNN_OUTPUTS_NUMBER(int _MNN_OUTPUTS_NUMBER){MNN_OUTPUTS_NUMBER = _MNN_OUTPUTS_NUMBER;}
	void setMNN_HIDDEN_UNITS(int _MNN_HIDDEN_UNITS){MNN_HIDDEN_UNITS = _MNN_HIDDEN_UNITS;}
	//print option file contents to standard output
	void write();
};
#endif
