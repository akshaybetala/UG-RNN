

#ifndef NN_h
#define NN_h 1
#include "Layer.h"
#include "General.h"
#include "NNperc.h"
#include <math.h>

// NN ver. 3.02 (27/11/2003)
// Copyright (C) Gianluca Pollastri 2003
//
// One-hidden layer Feedforward neural net.
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

class NN : public NNperc
{

  Layer_tanh* lower;


public:



  // Constructor. Parameters:
  // Number of input attributes, number of hidden units, number of output units;
  // t_NK contains the cardinalities of the attribute spans.

  NN(int t_NI, int t_NH, int t_NO, int* t_NK) :
    NNperc(t_NI, t_NH, t_NO, t_NK)
    {
      upper= new Layer_soft(NO,NK,0,NH);
      upper->set_output(1);
      upper->set_ninput(2);
      lower= new Layer_tanh(NH,NK,NI);
      lower->set_ninput(1);
    };

  // Constructor for a net with mixed inputs.
  // NI = number of input attributes (categorical inputs)
  // NIr = number of inputs (real valued)
  // ..
  // outp = output or non-output network (for backprop signal)
  // inp = input or non-inpput network (for backprop signal)


  NN(int t_NI,int t_NIr, int t_NH, int t_NO, int* t_NK,
	  int t_outp=1, int t_inp=1, int t_which=1) :
	NNperc(t_NI, t_NIr, t_NH, t_NO, t_NK, t_outp, t_inp)
    {
      upper= new Layer_soft(NO,NK,0,NH);

      if (outp)
        upper->set_output(1);
      upper->set_ninput(2);

      lower= new Layer_tanh(NH,NK,NI,NIr);
      if (inp)
        lower->set_ninput(1);
    };

  // Create/read a net from file
  NN(istream& is);
  void read(istream& is);

  // Forward pass
  void forward(int* I);
  void forward(Float* I);
  void forward(int* I1, Float* I2);
  void forward(Float* I1, Float* I2);


  // Backprop
  Float backward(Float* t, Float weight=1.0);
  Float* back_out() {return backprop;}

  // Update gradients
  void gradient(int* I, Float* t);
  void gradient(Float* I, Float* t);
  void gradient(int* I1, Float* I2, Float* t);
  void gradient(Float* I1, Float* I2, Float* t);
  
  // Update weights
  virtual void updateWeights(Float epsilon);
  virtual void updateWeightsL1(Float epsilon);
  virtual void updateWeightsClipped(Float epsilon);
  void resetGradient();
  virtual void initWeights(int seed);
  inline Float* out() { return upper->out(); };
  void write(ostream& os);

  void set_input(int vi) {
	  lower->set_ninput(vi);
	  inp=vi;
  }

  Float dlength() {
    return upper->dlength()+lower->dlength();
  }
};



#endif // NN_h
