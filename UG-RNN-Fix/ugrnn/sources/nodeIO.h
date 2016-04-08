#ifndef NODEIO_H
#define NODEIO_H
#include <string>
#include <vector>
#include <iostream>
#include "General.h"
using namespace std;
class NodeIO{
	protected:
	int categoricalInputSize;
	int realValuedInputSize;
	int outputSize;
	Float* I1;
	Float* I2;
	Float* out;

	public:
	//empty constructor
	//NodeIO();
	//constructor
	NodeIO(int _categoricalInputSize,int _realValuedInputSize,int _outputSize);
	//copy constructor
	NodeIO(NodeIO & _NodeIO);
	//assignment operator
	NodeIO operator=(NodeIO & _NodeIO);
	//destructor
	~NodeIO();
	//getters
	int getCategoricalInputSize(){return categoricalInputSize;}
	int getRealValuedInputSize(){return realValuedInputSize;}
	int getOutputSize(){return outputSize;}
	Float* getI1(){return I1;}
	Float* getI2(){return I2;}	
	Float getI1(int i){return I1[i];}
	Float getI2(int i){return I2[i];}
	Float getOut(int i){return out[i];}
	Float* getOut(){return out;}
	//setters
	//set element i to j
	void setI2(int i,Float j){I2[i]=j;}
	void setI2(){
		for(int i=0;i<realValuedInputSize;i++){
			I2[i]=0;	
		}
	}	
	void setOut(Float* _out){
		for(int i=0;i<outputSize;i++){
			out[i]=_out[i];
		}
	}
	//categorical init
	void categoricalInit(const string & symbol,const vector<string> & totalAtoms);
	//print I1
	void printI1(){
		for(int i=0;i<categoricalInputSize;i++){
			cout << I1[i] << " "; 	
		}
		cout << endl;
	}
	//print I2 
	void printI2(){
		for(int i=0;i<realValuedInputSize;i++){
			cout << I2[i] << " ";
		}
		cout << endl;
	}
	//print out
	void printOut(){
		for(int i=0;i<outputSize;i++){
			cout << out[i] << " ";
		}
		cout << endl;
	}
	
							
};
#endif
