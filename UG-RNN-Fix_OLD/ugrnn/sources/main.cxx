#include <math.h>
#include <string>
#include <vector>
#include "DataSet.h"
#include "myOptions.h"
#include "UGRNN.h"

using namespace std;

string floatToString(Float);

int main(int argc,char** argv){
//parsing input
string usageString("USAGE: ./ugrnnTrain -d dataset -t testID -f fold -m model -vl -tsval -restart epoch");
string datasetName("");
string testID("");
int fold=-1;
int model=-1;
int validate=0;
int testPlusVal=0;
int numberOfModels=10;
int restart = 0;
int epoch = 0;
if (argc<9) {
	cerr << usageString << endl;
	exit(1);
}
for(int i=0;i<argc;i++){
	if(strcmp(argv[i],"-d")==0){
		datasetName = argv[i+1];
	}
	else if(strcmp(argv[i],"-t")==0){
		testID = argv[i+1];
	} 
	else if(strcmp(argv[i],"-f")==0){
		fold = atoi(argv[i+1]);
	} 
	else if(strcmp(argv[i],"-m")==0){
		model = atoi(argv[i+1]);
	} 
	else if(strcmp(argv[i],"-vl")==0){		
		validate = 1;
	} 
	else if(strcmp(argv[i],"-tsval")==0){
		validate = 0;		
		testPlusVal = 1;
	}
	else if(strcmp(argv[i],"-restart")==0){
		restart = 1;
		epoch = atoi(argv[i+1]);
	} 
}
if(datasetName.compare("")==0 || testID.compare("")==0 || fold < 0 || model < 0 || epoch < 0 ){
	cerr << usageString;
	exit(1);
}
///////////////
//global variables
ifstream myOptStream;
stringstream mss;
stringstream fss;
myOptions myOpt;
string trainfile;
string datasetPath;
ifstream datasetStream;
ifstream trainStream;
ifstream validationStream;
ifstream testStream;
ofstream boundStream;
vector<int> trainVector;
vector<int> validationVector;
vector<int> testVector;
int trainLength=0;
int validationLength=0;
int testLength=0;
Molecule** train;
Molecule** validation;
Molecule** test;
vector<Float> bound(6,0);
Float CPU_time_usage;
//////////////////

//load myOptions file
mss << model;
myOptStream.open(("test/"+testID+"/fold0/model"+mss.str()+"/myOptions").c_str());
if(myOptStream.is_open()){
	myOpt = myOptions(myOptStream);	
}else{
	cerr << "CANNOT OPEN myOptions FILE" << endl;
	exit(1);
}
myOpt.write();
myOptStream.close();
/////////////////////
//load dataset
datasetPath = "./data/"+datasetName;
datasetStream.open(datasetPath.c_str());
if(!datasetStream.is_open()){
	cerr << "CANNOT OPEN DATASET " << datasetName << endl;	
	exit(1);
}
//loading dataset
clock_t begin = clock();
DataSet dataset(datasetStream,myOpt.getClasses());
datasetStream.close();
clock_t end = clock();
CPU_time_usage = get_CPU_time_usage(end,begin);
cout << "CPU_time_usage = " << CPU_time_usage << endl;
//cin >> dummy;
//////////////
//init UGRNN
fss << fold;
UGRNN* ugrnn = new UGRNN(dataset.mole,myOpt,dataset.length,"test/"+testID+"/fold"+fss.str()+"/model"+mss.str());
ugrnn->setFirstEpoch(0);
////////////
//init vector of molecules
trainStream.open(("test/"+testID+"/fold"+fss.str()+"/train").c_str());
if(!trainStream.is_open()){
	cerr << "CANNOT OPEN train FILE" << endl;
	exit(1);
}
trainStream >> trainLength;
for(int i=0;i<trainLength;i++){	
	int temp;
	trainStream >> temp;
	trainVector.push_back(temp);	
}
//debug
//vector<int>::iterator it = trainVector.begin();
//while(it!=trainVector.end()){
//	cout << *it++ << " ";
//}
//cout << endl;
/////// 
trainStream.close();
ugrnn->setTrainVector(trainVector);

if(testPlusVal == 1){
	testStream.open(("test/"+testID+"/fold"+fss.str()+"/test").c_str());
	if(!testStream.is_open()){
		cerr << "CANNOT OPEN test FILE" << endl;
   	}
	testStream >> testLength;
	for(int i=0;i<testLength;i++){	
		int temp;
		testStream >> temp;
		testVector.push_back(temp);	
	}
	//debug
	//vector<int>::iterator it = testVector.begin();
	//while(it!=testVector.end()){
	//	cout << *it++ << " ";
	//}
	//cout << endl;
	/////// 
	testStream.close();
	ugrnn->setTestVector(testVector);
	
	validationStream.open(("test/"+testID+"/fold"+fss.str()+"/validation").c_str());
	if(!validationStream.is_open()){
		cerr << "CANNOT OPEN validation FILE" << endl;
		exit(1);
	}
	validationStream >> validationLength;
	for(int i=0;i<validationLength;i++){	
		int temp;
		validationStream >> temp;
		validationVector.push_back(temp);	
	}

	//debug
	//vector<int>::iterator it = testVector.begin();
	//while(it!=testVector.end()){
	//	cout << *it++ << " ";
	//}
	//cout << endl;
	/////// 
	validationStream.close();
	ugrnn->setValidationVector(validationVector);
}

if(validate==1){
	validationStream.open(("test/"+testID+"/fold"+fss.str()+"/validation").c_str());
	if(!validationStream.is_open()){
		cerr << "CANNOT OPEN validation FILE" << endl;
		exit(1);
	}
	validationStream >> validationLength;
	for(int i=0;i<validationLength;i++){	
		int temp;
		validationStream >> temp;
		validationVector.push_back(temp);	
	}

	//debug
	//vector<int>::iterator it = testVector.begin();
	//while(it!=testVector.end()){
	//	cout << *it++ << " ";
	//}
	//cout << endl;
	/////// 
	validationStream.close();
	ugrnn->setValidationVector(validationVector);
}

////////////////////////////////////
//init train validate and test 
train = new Molecule*[trainLength];
test = new Molecule*[testLength];
validation = new Molecule*[validationLength];
for(int i=0;i<trainLength;i++){
	train[i] = dataset.mole[trainVector[i]];
}
if(validate == 1){
	for(int i=0;i<validationLength;i++){
		validation[i] = dataset.mole[validationVector[i]];
	}
}
if(testPlusVal == 1){
	for(int i=0;i<validationLength;i++){
		validation[i] = dataset.mole[validationVector[i]];
	}
	for(int i=0;i<testLength;i++){
		test[i] = dataset.mole[testVector[i]];
		//debug
		//cout << "test target " << i << " " << test[i]->target << endl;
		///////
	}
}
/////////////////////

//normalize targets

dataset.normalizeTargets(train,trainLength);
bound[0] = dataset.min;
bound[1] = dataset.max;
bound[2] = dataset.myMin;
bound[3] = dataset.myMax;
bound[4] = dataset.minLogP;
bound[5] = dataset.maxLogP;

boundStream.open(("test/"+testID+"/fold"+fss.str()+"/model"+mss.str()+"/bound").c_str());
vector<Float>::iterator it = bound.begin();
boundStream << floatToString(6) << endl;
while(it != bound.end()){
	boundStream << floatToString(*it++) << endl;	
}
boundStream.close();
//preprocessing ends
 
//restart
if(restart == 1){
	string modelDir = "test/"+testID+"/fold0/model"+mss.str();
	stringstream es;
	es << epoch;	
	ugrnn->read(modelDir+"/ONN"+es.str(),modelDir+"/MNN"+es.str());
	ugrnn->setFirstEpoch(epoch);
}
/////////
//run train, validate and test

if(validate == 1){
	dataset.normalizeTargets(validation,validationLength,dataset.max,dataset.min);
	ugrnn->trainAndValidate(train,validation,bound,numberOfModels);	
}
else if(testPlusVal == 1){
	dataset.normalizeTargets(test,testLength,dataset.max,dataset.min);
	dataset.normalizeTargets(validation,validationLength,dataset.max,dataset.min);	
	ugrnn->trainAndValidate(train,validation,bound,numberOfModels);	
	ugrnn->predict(test,bound,numberOfModels);	
	//ugrnn->printDescriptors(train,test,bound,numberOfModels);
}
else{
	cerr << usageString << endl;
    exit(0);	
}
///////////////////
//free memory
//unload graph
//for(int i=0;i<dataset.length;i++){
//	dataset.mole[i]->g->unloadGraph();
//}
//////////////
delete ugrnn;
delete[] train;
delete[] test;
delete[] validation;
/////////////

}

string floatToString(Float f){
	stringstream ss;
	ss << f;
	return ss.str();
}
