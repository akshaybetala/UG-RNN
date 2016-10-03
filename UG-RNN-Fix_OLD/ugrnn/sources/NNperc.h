

#ifndef NNperc_h
#define NNperc_h 1
#include "Layer.h"
#include "General.h"
#include <math.h>

// NN ver. 3.02 (27/11/2003)
// Copyright (C) Gianluca Pollastri 2003
//
// one layer
// Input: categorical data (one-hot), real valued or mixed.
// ouput: softmax.
// Cost function: the proper one
// Gradient: plain backpropagation (no momentum)
//
// 
// Version 3.0:
// compatible with Layer 3.0
//
// Version 3.02
// Added updateWeightsL1

class NNperc
{

protected:
  int NI;
  int NIr;
  int NItot;
  int NH;
  int NO;
  int* NK;
  int* NK2;
  int which;
  int outp;
  int inp;

  Float* backprop;

  Layer_soft* upper;

public:

  // Constructor. Parameters:
  // Number of input attributes, number of hidden units, number of output units;
  // t_NK contains the cardinalities of the attribute spans.

  NNperc(int t_NI, int t_NH, int t_NO, int* t_NK) :
    NI(t_NI), NH(t_NH), NO(t_NO), NK(t_NK)
    {
      NK=new int[NI];
      NItot=0;
      for (int i=0; i<NI; i++) {
	NK[i]=t_NK[i];
	NItot += NK[i];
      }
      upper= new Layer_soft(NO,NK,NI);
      upper->set_output(1);
      upper->set_ninput(1);

      NIr=0;
      outp=1;
      inp=1;
    };

  // Constructor for a net with mixed inputs.
  // NI = number of input attributes (categorical inputs)
  // NIr = number of inputs (real valued)
  // ..
  // outp = output or non-output network (for backprop signal)
  // inp = input or non-inpput network (for backprop signal)


  NNperc(int t_NI,int t_NIr, int t_NH, int t_NO, int* t_NK,
	  int t_outp=1, int t_inp=1, int t_which=1) :
	NI(t_NI), NIr(t_NIr), NH(t_NH), NO(t_NO), NK(t_NK), outp(t_outp), inp(t_inp)
    {
      int i;
      NK=new int[NI];
      NItot=0;
      for (i=0; i<NI; i++) {
	NK[i]=t_NK[i];
	NItot += NK[i];
      }
      NK2=new int[NIr];
      for (i=0; i<NIr; i++)
	NK2[i]=1;

      which=1;
      upper= new Layer_soft(NO,NK,NI);

      if (outp)
        upper->set_output(1);
      upper->set_ninput(1);

      backprop=new Float[NItot+NIr]; 	
    };

  // Create/read a net from file
  NNperc(istream& is);
  virtual void read(istream& is);

  // Forward pass
  virtual void forward(int* I);
  virtual void forward(Float* I);
  virtual void forward(int* I1, Float* I2);
  virtual void forward(Float* I1, Float* I2);


  Float f_cost(Float* t) {
	return upper->f_cost(t);
  }

  // Backprop
  virtual Float backward(Float* t, Float weight=1.0);
  Float* back_out() {return backprop;}

  // Update gradients
  virtual void gradient(int* I, Float* t);
  virtual void gradient(Float* I, Float* t);
  virtual void gradient(int* I1, Float* I2, Float* t);
  virtual void gradient(Float* I1, Float* I2, Float* t);
  
  // Update weights
  virtual void updateWeights(Float epsilon);
  virtual void updateWeightsL1(Float epsilon);
  virtual void updateWeightsClipped(Float epsilon);
  virtual void resetGradient();
  virtual void initWeights(int seed);
  inline Float* out() { return upper->out(); };
  virtual void write(ostream& os);

  virtual void set_input(int vi) {
	  upper->set_ninput(vi);
	  inp=vi;
  }
  void set_output(int vo) {
	  upper->set_output(vo);
	  outp=vo;
  }


  inline int get_NI() { return NI; };
  inline int get_NIr() { return NIr; };
  inline int get_NO() { return NO; };
  inline int get_NH() { return NH; };

  Float dlength() {
    return upper->dlength()+upper->dlength();
  }
};



#endif // NN_h
