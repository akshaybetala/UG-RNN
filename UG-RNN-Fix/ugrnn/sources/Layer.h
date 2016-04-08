
#ifndef Layer_h
#define Layer_h 1
#include "General.h"
#include <math.h>



// Layer ver 3.03
// 12/12/2003
// Copyright (C) Gianluca Pollastri 2003
//
// ANN Layers
// Linear, tanh and softmax outputs
// Categorical (one-hot), real-valued and mixed inputs.
//
//
// In version 3.0:
// -fixed 'overflow problem' that output error=0
//  for saturated softmax units
// -all (but the last one..) compatibility issues fixed
// -added updateWeightsClipped
//
// In version 3.01
// -fixed all the versions of gradient();
// -NUtot made an attribute
//
// In version 3.02:
// -fixed Linux bug in initWeights
// -added updateWeightsL1
//
// In version 3.03:
// -gradient fixed for softmax: (y-t)x, no f1


class Layer
{
protected:
int NY;
int NU;
int NUr;
int* NK;

Float* Y;


Float* A;
Float* U;	 //NU*NK
Float* delta;    //NY
Float* backprop; //NU*NK

Float*** W;
Float*** dW;

Float* B;        //NY
Float* dB;       //NY

int output; 	//0=no,1=yes
int ninput;		//0=input layer,1=just real side backprop,2=full backprop

int NUtot,NUplain;


void alloc(int NY, int NU, int* NK);


public:

void softmax();
void squash();

// Destructor
virtual ~Layer(){
  delete[] NK;
  delete[] Y;
  delete[] A;
  delete[] U;
  delete[] delta;
  delete[] backprop;
  delete[] B;
  delete[] dB;
  for (int y=0; y<NY; y++) {
    for (int u=0; u<NU+NUr; u++) {
      delete[] W[y][u];
      delete[] dW[y][u];
    }
    delete[] W[y];
    delete[] dW[y];
  }
  delete[] W;
  delete[] dW; 
}
// Constructor
// Categorical inputs

Layer(int t_NY, int* t_NK, int t_NU) :
	NY(t_NY), NU(t_NU)
{
NUr=0;
NK=new int[NU];
for (int i=0; i<NU; i++)
	NK[i]=t_NK[i];
alloc(NY,NU,NK);
ninput=0;
output=0;

NUplain =0;
for (int u=0; u<NU; u++) {
  NUplain += NK[u];
}
}

// Constructor
// Real-valued inputs

Layer(int t_NY, int t_NU) :
	NY(t_NY), NU(t_NU)
{
NK=new int[NU];
for (int i=0; i<NU; i++)
	NK[i]=1;
NUr=0;
alloc(NY,NU,NK);
ninput=0;
output=0;

NUplain =0;
for (int u=0; u<NU; u++) {
  NUplain += NK[u];
}
//cout << NUplain << " " << NUtot << "a " << flush;
}

// Constructor
// Mixed inputs (NU categorical attributes, NUr real-valued)

Layer(int t_NY, int* t_NK, int t_NU, int t_NUr) :
	NY(t_NY), NU(t_NU), NUr(t_NUr)
{
int i;
NK=new int[NU+NUr];
for (i=0; i<NU; i++)
	NK[i]=t_NK[i];
for (i=NU; i<NU+NUr; i++)
	NK[i]=1;
alloc(NY,NU+NUr,NK);
ninput=0;
output=0;


NUplain =0;
for (int u=0; u<NU; u++) {
  NUplain += NK[u];
}
//cout << NY << " " << flush;
//cout << NUplain << " " << NUtot << "b " << flush;
}

Layer(istream& is);



void set_ninput(int vi) {
  ninput=vi;
};

void set_output(int vo) {
  output=vo;
};


void read(istream& is);
void write(ostream& os);

virtual void forward(int* I);
virtual void forward(Float* I);
virtual void forward(int* I1, Float* I2);
virtual void forward(Float* I1, Float* I2);

virtual Float f1(int y);
virtual Float f_cost(Float* t);
Float log_cost(Float* t);
Float sq_cost(Float* t);

virtual Float backward(Float* t, Float weight=1.0);

void gradient(int* I);
void gradient(Float* I);
void gradient(int* I1, Float* I2);
void gradient(Float* I1, Float* I2);
void scaleGradient(int scaleFactor);
void gradient();

void getDWcontribution();
void getDWcontribution(Float* I1, Float* I2);

virtual void updateWeights(Float epsilon);
virtual void updateWeights(Float epsilon, Float factor);
virtual void updateWeightsL1(Float epsilon);
virtual void updateWeightsClipped(Float epsilon);
virtual void updateWeightsClippedWeighted(Float epsilon,Float lambda);
void resetGradient();
virtual void initWeights(int seed);

inline Float* back_out() { return backprop; }
inline Float* Aout() { return A; }
inline Float* out() { return Y; }


inline int get_NY() { return NY; }
inline int get_NU() { return NU; }
inline int* get_NK() { return NK; }

inline Float*** get_dW() { return dW; }

Float dlength();

void set_dW(Float*** newdW);


};


class Layer_tanh : public Layer
{

public:


Layer_tanh(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
}

Layer_tanh(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
}

Layer_tanh(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
}

Layer_tanh(istream& is) :
Layer(is)
{
}

void forward(int* I)
{
Layer::forward(I);
squash();
}

void forward(Float* I)
{
Layer::forward(I);
squash();
}

void forward(int* I1,Float* I2)
{
Layer::forward(I1,I2);
squash();
}

void forward(Float* I1,Float* I2)
{
Layer::forward(I1,I2);
squash();
}


Float backward(Float* t, Float weight=1.0)
{
return Layer::backward(t,weight);
}


Float f1(int y);

Float f_cost(Float* t);


void initWeights(int seed)
{
Layer::initWeights(seed);
}

void updateWeights(Float epsilon)
{
Layer::updateWeights(epsilon);
}

void updateWeights(Float epsilon, Float factor)
{
Layer::updateWeights(epsilon, factor);
}

};



class Layer_soft : public Layer
{

public :

Layer_soft(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
}

Layer_soft(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
}

Layer_soft(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
}

Layer_soft(istream& is) :
Layer(is)
{
}

void forward(int* I)
{
Layer::forward(I);
softmax();
}

void forward(Float* I)
{
Layer::forward(I);
softmax();
}

void forward(int* I1,Float* I2)
{
Layer::forward(I1,I2);
softmax();
}

void forward(Float* I1,Float* I2)
{
Layer::forward(I1,I2);
softmax();
}


Float backward(Float* t, Float weight=1.0);

Float f1(int y);

Float f_cost(Float* t);

void initWeights(int seed)
{
Layer::initWeights(seed);
}

void updateWeights(Float epsilon)
{
Layer::updateWeights(epsilon);
}

};




/*

class Selector : public Layer
{
int* edges;

public:


Selector(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
edges=new int[NY];
}

Selector(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
edges=new int[NY];
}

Selector(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
edges=new int[NY];
}

Selector(istream& is) :
Layer(is)
{
edges=new int[NY];
}


void forward(int* I)
{}
void forward(Float* I);
void forward(int* I1, Float* I2)
{}
void forward(Float* I1, Float* I2)
{}

Float f1(Float a);

Float f_cost(Float* t);

Float backward(Float* t, Float weight=1.0);

void setEdges();

void initWeights(int seed)
{
Layer::initWeights(seed);
setEdges();
}

void updateWeights(Float epsilon)
{
Layer::updateWeights(epsilon);
setEdges();
}

};



*/

#endif
