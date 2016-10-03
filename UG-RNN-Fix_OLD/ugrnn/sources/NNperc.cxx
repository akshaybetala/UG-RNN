
// NN ver. 3.02 (27/11/2003)
//
// Copyright (C) Gianluca Pollastri 2003


#include "NNperc.h"


NNperc::NNperc(istream& is)
{
  is >> NO >> NH >> NI >> NIr >> which >> outp >> inp;
  upper= new Layer_soft(is);
  if (outp)
    upper->set_output(1);
  upper->set_ninput(inp);

  int i;
  NK=new int[NI];
  NItot=0;
  for (i=0; i<NI; i++) {
	NK[i]=upper->get_NK()[i];
	NItot += NK[i];
  }
  NK2=new int[NIr];
  for (i=0; i<NIr; i++)
	NK2[i]=upper->get_NK()[NI+i];
  backprop=new Float[NItot+NIr];
}


void
NNperc::read(istream& is)
{
  is >> NO >> NH >> NI >> NIr >> which >> outp >> inp;
  upper->read(is);
  if (outp)
  	upper->set_output(1);
  upper->set_ninput(inp);

  int i;
  NItot =0;
  for (i=0; i<NI; i++) {
	NK[i]=upper->get_NK()[i];
	NItot += NK[i];
  }
  for (i=0; i<NIr; i++)
	NK2[i]=upper->get_NK()[NI+i];
}



void
NNperc::forward(int* I)
{
  upper->forward(I);
}

void
NNperc::forward(Float* I)
{
  upper->forward(I);
}

void
NNperc::forward(int* I1,Float* I2)
{
  upper->forward(I1,I2);
}

void
NNperc::forward(Float* I1,Float* I2)
{
  upper->forward(I1,I2);
}

Float
NNperc::backward(Float* t, Float weight)
{
  Float err=upper->backward(t,weight);
  if (inp==1)
    for (int r=NItot;r<NItot+NIr;r++)
      backprop[r]=upper->back_out()[r];
  else if (inp==2)
    for (int r=0;r<NItot+NIr;r++)
      backprop[r]=upper->back_out()[r];
  return err;
}

void
NNperc::gradient(int* I, Float* t)
{
  upper->gradient(I);
}
void
NNperc::gradient(Float* I, Float* t)
{
  upper->gradient(I);
}
void
NNperc::gradient(int* I1,Float* I2, Float* t)
{
  upper->gradient(I1,I2);
}
void
NNperc::gradient(Float* I1,Float* I2, Float* t)
{
  upper->gradient(I1,I2);
}



void
NNperc::resetGradient()
{
  upper->resetGradient();
}

void
NNperc::updateWeights(Float epsilon)
{
  upper->updateWeights(epsilon);
}

void
NNperc::updateWeightsL1(Float epsilon)
{
  upper->updateWeightsL1(epsilon);
}

void
NNperc::updateWeightsClipped(Float epsilon)
{
  upper->updateWeightsClipped(epsilon);
}

void
NNperc::initWeights(int seed)
{
  upper->initWeights(seed);
}


void
NNperc::write(ostream& os)
{
  os << NO << " " << NH<< " " << NI<< " " << NIr <<" ";
  os << which << " " << outp << " " << inp << "\n";
  upper->write(os);
}



