#include <math.h>
#include <string>
#include <vector>
#include "DataSet.h"
#include "myOptions.h"
#include "UGRNN.h"

using namespace std;

void parseInput(char,char**,string &,string &,string &,string &);
vector<Float> initBound(ifstream &);
template <typename T>
void printVector(const vector<T> &);


int main(int argc,char** argv){
//variables
string inputString("");
string modelDir("");
string fold("0");
string model("0");
UGRNN* ugrnn;
myOptions myOpt;
vector<Float> bound;
vector<Float> predictedValues;

//parse input
parseInput(argc,argv,inputString,modelDir,fold,model);

//read myOptions
ifstream myOptStream;
string str = modelDir+"/fold0/model"+model+"/myOptions";
myOptStream.open((str).c_str());
if(myOptStream.is_open()){
	myOpt = myOptions(myOptStream);	
}else{
	cerr << "CANNOT OPEN myOptions FILE\n" << str << endl;
	exit(1);
}
myOpt.write();
myOptStream.close();

//load input file
ifstream inputStream;
inputStream.open((inputString).c_str());
if(!inputStream.is_open()){
	cerr << "CANNOT OPEN DATASET " << inputString << endl;	
	exit(1);
}
DataSet input(inputStream,myOpt.getClasses());
inputStream.close(); 

//init UGRNN
ugrnn = new UGRNN(input.mole,myOpt,input.length,modelDir+"/fold"+fold+"/model"+model,1);
vector<int> inputVector;
for(int i=0;i<input.length;i++){
	inputVector.push_back(i);
}
ugrnn->setInputVector(inputVector);

//predict
ifstream boundStream;
boundStream.open((modelDir+"/fold0/model"+model+"/bound").c_str());
if(!boundStream.is_open()){
	cerr << "CANNOT OPEN BOUND FILE" << endl;
	exit(1);	
}
bound = initBound(boundStream);
//debug
//printVector(bound);
///////
predictedValues = ugrnn->predictTarget(input.mole,bound,10);
//debug
printVector(predictedValues);
///////

//free memory
//unload graph
//for(int i=0;i<input.length;i++){
//	input.mole[i]->g->unloadGraph();
//}
//////////////
delete ugrnn;

}




void parseInput(char argc,char** argv,string & inputString,string & modelDir,string & fold,string & model){
	string usageString("USAGE: ./ugrnnPredict -i input -md modelDirectory -f fold -m model\n");	
	for(int i=0;i<argc;i++){
		if(strcmp(argv[i],"-i")==0){
			inputString = argv[i+1];
		}
		else if(strcmp(argv[i],"-md")==0){
			modelDir = argv[i+1];
		}
		else if(strcmp(argv[i],"-f")==0){
			 fold = argv[i+1];
		}
		else if(strcmp(argv[i],"-m")==0){
			model = argv[i+1];
		} 
	}
	if(inputString.compare("")==0 || modelDir.compare("")==0){
		cerr << usageString;
		exit(1);
	}	
}

template <typename T>
void printVector(const vector<T> & vect){
	typename vector<T>::const_iterator it = vect.begin();
	while(it!=vect.end()){
		cout << *it++ << " ";	
	}
	cout << endl;	
}

vector<Float> initBound(ifstream & boundStream){
		vector<Float> out;
		int numberOfBounds; 
		int counts = 0;
		Float bound;		
		
		boundStream >> numberOfBounds; 		
		while(counts < numberOfBounds){
			boundStream >> bound;
			out.push_back(bound);
			boundStream;
			counts++;	
		}
		return out;
}

