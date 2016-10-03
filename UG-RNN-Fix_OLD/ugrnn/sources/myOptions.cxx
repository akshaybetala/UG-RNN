#include"myOptions.h"
myOptions::myOptions(ifstream& is){
	char* str = new char[1024];
	while(is){
		is >> str;
		if(!strcmp(str,"classes")) is >> classes;
		else if(!strcmp(str,"ONN_OUTPUTS_NUMBER")) is >> ONN_OUTPUTS_NUMBER;
		else if(!strcmp(str,"ONN_HIDDEN_UNITS")) is >> ONN_HIDDEN_UNITS;
		else if(!strcmp(str,"MNN_OUTPUTS_NUMBER")) is >> MNN_OUTPUTS_NUMBER;
		else if(!strcmp(str,"MNN_HIDDEN_UNITS")) is >> MNN_HIDDEN_UNITS;
		else if(!strcmp(str,"seed")) is >> seed;
		else if(!strcmp(str,"epsilon")) is >> epsilon;
		else if(!strcmp(str,"nEpochs")) is >> nEpochs;
		else if(!strcmp(str,"batchBlocks")) is >> batchBlocks;
		else if(!strcmp(str,"folds")) is >> folds;
		else if(!strcmp(str,"epsilonMin")) is >> epsilonMin;
		else if(!strcmp(str,"epsEpoch")) is >> epsEpoch;
		else if(!strcmp(str,"properties")) is >> properties;
	}
	if(classes == 0){
		cerr << "Classes = 0" << " program will quit" << endl;
		exit(1);
	}
	delete [] str;	
	
};

myOptions::~myOptions(){}
void myOptions::write(){
	cout << "///////////////////////////////////////////////////////" << endl;
        cout << "OPTIONS" << endl;
	cout << "classes " << getClasses() << endl;
	cout << "ONN_OUTPUTS_NUMBER " << getONN_OUTPUTS_NUMBER() << endl;
	cout << "ONN_HIDDEN_UNITS " << getONN_HIDDEN_UNITS() << endl;
	cout << "MNN_OUTPUTS_NUMBER " << getMNN_OUTPUTS_NUMBER() << endl;
	cout << "MNN_HIDDEN_UNITS " << getMNN_HIDDEN_UNITS() << endl;
	cout << "seed " << getSeed() << endl;
	cout << "epsilon " << getEpsilon() << endl;
	cout << "nEpochs " << getNEpochs() << endl;
	cout << "batchBlocks " << getBatchBlocks() << endl;
	cout << "folds " << getFolds() << endl;
	cout << "epsilonMin " << getEpsilonMin() << endl;
	cout << "epsEpoch " << getEpsEpoch() << endl;
	cout << "properties " << getProperties() << endl;
	cout << "///////////////////////////////////////////////////////" << endl;
}
