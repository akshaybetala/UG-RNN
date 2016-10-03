#include "nodeIO.h"
NodeIO::NodeIO(int _categoricalInputSize,int _realValuedInputSize,int _outputSize):categoricalInputSize(_categoricalInputSize),realValuedInputSize(_realValuedInputSize),outputSize(_outputSize){
	//allocating memory space
	I1 = new Float[_categoricalInputSize];
	memset(I1,0,_categoricalInputSize*sizeof(Float));
	I2 = new Float[_realValuedInputSize];
	memset(I2,0,_realValuedInputSize*sizeof(Float));
	out = new Float[_outputSize];
	memset(out,0,_outputSize*sizeof(Float));
		
}

NodeIO::NodeIO(NodeIO& _NodeIO){
	cout << "calling NodeIO copy constructor" << endl;
}
NodeIO::~NodeIO(){
	//cout << "calling NoedIO destructor" << endl;
	delete[] I1;
	delete[] I2;
	delete[] out;
}

void NodeIO::categoricalInit(const string & symbol,const vector<string> & totalAtoms){
	for(int i=0;i<totalAtoms.size();i++){
		if(symbol.compare(totalAtoms[i])==0){
                	I1[i]=1.0;
                        break;
                }
        }
}


