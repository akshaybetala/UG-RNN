//////////////////////////////////////////////////////////////////////////////////////
// Factor.h  --  class definition for matlab-compatible factor class
//
// A few functions are defined only for MEX calls (construction & load from matlab)
// Most others can be used more generally.
// 
//////////////////////////////////////////////////////////////////////////////////////
//
// Written by Alex Ihler
// Copyright (C) 2010 Alexander Ihler; distributable under GPL -- see README.txt
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef __FACTOR_H
#define __FACTOR_H

#include "Variables.h"
#include "subindex.h"

//#include <assert>
#include <cmath>
#include <iostream>
#include <limits>
#include <float.h>

#include <vector>

class Factor {
 public:

  /////////////////////////////
  // Constructors / Destructors
  /////////////////////////////

  Factor() : v_() { mxInit(); N_=1; t_=new value[1]; t_[0]=1.0; };				// blank constructor
  Factor(Factor const& copy) : v_() { t_=NULL; N_=0; mxInit(); *this = copy; };	// copy constructor
  Factor(const value s) : v_() { mxInit(); N_=1; t_=new value[1]; t_[0]=s; };	// scalar constructor
  Factor(Variables const& Vars, const value Scalar=0.0) : v_(Vars) { 			// scalar "sized" const.
		mxInit(); N_=calcNumEl(); t_=new value[N_]; 
		for (vsize i=0;i<N_;++i) t_[i]=Scalar; 
  };

  ~Factor() { 
	if (!mxAvail()) {					// if we exist as a matlab class, do nothing (?)
	  if (t_) delete[] t_;				// if not, deallocate our memory
	}
  }

  // Assignments & copy constructors
  Factor& operator=(Factor const& rhs) {	
    if (this!=&rhs) {
	  if (mxAvail()) {
#ifdef MEX		// Never executed (mxAvail false) if non-MEX code ; functions don't exist
	    if (v_ == rhs.v_) for (vsize i=0;i<N_;i++) t_[i]=rhs.t_[i];
		else {
		  mxArray* vNew = mxCreateNumericMatrix(1, rhs.nvar(), mxUINT32_CLASS,mxREAL);
 		  mxDestroyArray(mxGetField(mHandle,0,"t"));						// release old table & make new
		  mxArray* tNew = mxCreateNumericArray(rhs.nvar(),rhs.dims(), mxDOUBLE_CLASS,mxREAL);
		  t_=(value*)mxGetData(tNew); 
		  v_.mxSet(vNew,mxGetDimensions(tNew)); v_=rhs.v_; 					// copy off variable data
		  N_=calcNumEl(); for (vsize j=0;j<N_;j++) t_[j]=rhs.t_[j];			// copy off table elements
		  mxSetField(mHandle,0,"v",vNew); mxSetField(mHandle,0,"t",tNew);	// link matrices into structure
		}	
#endif
	  } else {
		 v_ = rhs.v_; 														// copy off variables
	     if (t_!=NULL) delete[] t_;											// release old table of values
	     N_=rhs.numel();  t_=new value[N_]; 								// allocate for new table
		 for (vsize j=0;j<N_;j++)  t_[j]=rhs.t_[j];							// and copy over
	  }
	}
	return *this;
  }

#ifdef MEX
/////////////////////////////////////////////////////////////////////////////////////////
// MEX Class Wrapper Functions
//   mxAvail()        : true if we already represent a matlab object
//   mxCopy(mxArray*) : create a *copy* of the matlab data & wrap it
//   mxSet(mxArray*)  : wrap passed matlab data; assumes irrevocable, exclusive ownership 
//   mxNew()          : create a brand-new matlab class object and copy current data in
//   mxGet()          : Return a pointer to the wrapped object (created if required)
/////////////////////////////////////////////////////////////////////////////////////////

  // mxCopy : Create a *copy* of the matlab data and wrap it
  Factor& mxCopy(mxArray const* mNew) {	
	//std::cout<< "Copying factor for wrap...\n";
	mxArray* mCopy = mxDuplicateArray(mNew);			// duplicate the data (deep copy)
	mxSet(mCopy);										// and wrap the copy
  }
  
  // mxSet : Wrap a given matlab memory structure
  Factor& mxSet(mxArray* mNew) {						// make object wrap a matlab class object
	//std::cout<<"Setting factor wrapper ("<<mxAvail()<<")\n";
	if (mxAvail()) mxDestroyArray(mHandle);				//   destroy old matlab object if it existed
	mHandle=mNew;
    mxArray* mxt=mxGetField(mNew, 0, "t");				// wrap the table matrix
	if (mxt==NULL) mexErrMsgTxt("Error converting factor argument");
    N_ = mxGetNumberOfElements(mxt);
    t_ = (value*) mxGetData(mxt);
  	mxArray* mxv=mxGetField(mNew, 0, "v");				// ask v_ to wrap the variable matrix
	if (mxv==NULL) mexErrMsgTxt("Error converting factor argument");
    v_.mxSet(mxv, mxGetDimensions(mxt)); 
  }

  // mxNew : Create a new matlab memory structure to hold our current data
  Factor& mxNew() {
	mxArray* mNew;
	// Create new empty matlab factor from blank constructor
    int retval = mexCallMATLAB(1,&mNew,0,NULL,"factor");
	if (retval) mexErrMsgTxt("Error creating new factor");
	// If we currently have a non-empty factor, we need to create things manually
	if (nvar()>0) {
		mHandle=mNew;
	  // Create table storage and copy in table data
      mxArray* tNew = mxCreateNumericArray(nvar(),dims(), mxDOUBLE_CLASS,mxREAL);
      value* tNew_ = (value*)mxGetData(tNew); for (vsize i=0;i<N_;i++) tNew_[i]=t_[i];
	  if (mxAvail()) mxDestroyArray(mxGetField(mHandle,0,"t")); else delete[] t_;
	  mxSetField(mHandle,0,"t",tNew);
	  // Ditto for variables
	  v_.mxNew(mxGetDimensions(tNew)); 
	  mxSetField(mHandle,0,"v",v_.mxGet());
  } else {	// if we have a scalar factor, just copy in the correct value
		value v=t_[0]; mxSet(mNew); t_[0]=v;
	}
  };
	
  // mxAvail : Check whether we are currently a matlab matrix
  inline bool mxAvail() const { return mHandle!=NULL; };

  // mxGet   : Return current matlab handle (create if doesn't exist)
  mxArray* mxGet(void) { 
	//std::cout<<"Getting wrapped factor ("<<mxAvail()<<")\n";
	if (!mxAvail()) mxNew();
	return mHandle; 
  };
#else
  inline bool mxAvail() const { return false; };
#endif
////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////
  // Accessor Functions  
  /////////////////////////////
  const vindex     nvar()      const { return v_.nvar();  };// # of variables
  const Variables& vars()      const { return v_;   };		// variable IDs
  const Variables& variables() const { return v_;   };		// (same as "vars")
  const vsize*     dims()      const { return v_.dims();};  // variable dimensions
  const value*     table()     const { return t_;   };		// table of factor values

  //////////////////////////////////////
  // Boolean checks on factor properties
  //////////////////////////////////////
  bool isempty()  const { return (t_==NULL); };			// empty factor
  bool isnan()    const { bool b=false; for (vsize i=0;i<N_;i++) b|=isnan(t_[i]);    return b; };
  bool isfinite() const { bool b=true;  for (vsize i=0;i<N_;i++) b&=isfinite(t_[i]); return b; };
  bool isscalar() const { return (N_==1); };
  vsize numel()   const { return N_; };
  vsize calcNumEl() const { vsize n=1; vsize const *d; d=dims(); for (vindex i=0;i<nvar();i++) n*=d[i]; return (n>1)?n:1; };
  //vsize numel()   const { vsize n=1; vindex N=v.nvar(); for (vindex i=0;i<N;i++) n*=dims[i]; return n; };

  //////////////////////////////////////
  // Direct value accessor
  //////////////////////////////////////
  value operator[] (vsize v) const { return t_[v]; };	// access table elements
  value& operator[] (vsize v)      { return t_[v]; };	// rewrite table elements

////////////////////////////////////////////////////////////////////////////////
// Unary transformations (in-place); see outside class def'n for const versions
////////////////////////////////////////////////////////////////////////////////
  inline Factor& abs(void)           { unaryOp( unOpAbs()   ); return *this; };
  inline Factor& exp(void)           { unaryOp( unOpExp()   ); return *this; };
  inline Factor& log(void)           { unaryOp( unOpLog()   ); return *this; };
  inline Factor& log2(void)          { unaryOp( unOpLogL(std::log(2)) ); return *this; };
  inline Factor& log10(void)         { unaryOp( unOpLog10() ); return *this; };
//  inline Factor& power(double pow)   { unaryOp( unOpPower(pow) ); return *this; };

  // Functors defined for unary operations (transformations) to the factor table
//  struct unOpPower { value b; unOpPower(value bb):b(bb){};  value operator()(value a) { return std::pow(a,b); } };
  struct unOpAbs   { value operator()(value a)  { return std::abs(a);  } };
  struct unOpExp   { value operator()(value a)  { return std::exp(a);  } };
  struct unOpLog   { value operator()(value a)  { return std::log(a);  } };
  struct unOpLogL  { value l; unOpLogL(value L):l(L) {}; value operator()(value a)  { return std::log(a)/l;  } };
  struct unOpLog10 { value operator()(value a)  { return std::log10(a);} };
  template<typename Function> void unaryOp( Function Op ) {
	for (vsize i=0;i<N_;i++) t_[i] = Op(t_[i]);
  }

////////////////////////////////////////////////////////////////////////////////
// Basic factor operations (+,-,*,/)
////////////////////////////////////////////////////////////////////////////////

  // Binary operations (eg A + B); returns new object
  template<typename Function> Factor binaryOp( const Factor B, Function Op) const {
	Variables v = v_ + B.v_;  						// expand scope to union
	Factor F(v);             						//  and create target factor
	subindex s1(v,v_), s2(v,B.v_); 					// index over A and B & do the op
	for (vindex i=0; i<F.N_; ++i,++s1,++s2) F[i]=Op(t_[s1], B[s2]);
	return F; 										// return the new copy
  }
  template<typename Function> Factor binaryOp( const value B, Function Op) const {
	Factor F=*this; F.binaryOpIP(B,Op); return F; // for scalar args, define with an in-place operator
  }
  // Binary in-place operations (eg A += B); returns reference to modified A
  template<typename Function> Factor& binaryOpIP( const Factor B, Function Op) {
	Variables v = v_ + B.v_;  								// expand scope to union
	if (v != v_) *this = binaryOp(B,Op); 					// if A's scope is too small, call binary op
	else {                         
	  subindex s2(v_,B.v_);       								// otherwise create index over B
	  for (vindex i=0; i<N_; ++i,++s2) Op.IP(t_[i] , B[s2]);	// and do the operations
	}                                               
	return *this; 
  }
  template<typename Function> Factor& binaryOpIP( const value B, Function Op) {
	for (vsize i=0;i<N_;i++) Op.IP(t_[i] , B); return *this;	// simplifies for scalar args
  }

  // Functors defined for binary operations on the factor table : Op(a,b) and Op.IP(a,b) (in-place version)
  struct binOpPlus   { 
	  value  operator()(value  a, const value b) { return a+b; }; 
	  value&         IP(value& a, const value b) { return a+=b;}; 
  };
  struct binOpMinus  { 
	  value  operator()(value  a, const value b) { return a-b; }; 
	  value&         IP(value& a, const value b) { return a-=b;}; 
  };
  struct binOpTimes  { 
	  value  operator()(value  a, const value b) { return a*b; }; 
	  value&         IP(value& a, const value b) { return a*=b;}; 
  };
  struct binOpDivide { 
	  value  operator()(value  a, const value b) { return (b) ? a/b  : 0;  }; 
	  value&         IP(value& a, const value b) { return (b) ? a/=b : a=0;}; 
  };
  struct binOpPower { 
	  value  operator()(value  a, const value b) { return std::pow(a,b);  }; 
	  value&         IP(value& a, const value b) { return a=std::pow(a,b);}; 
  };

  Factor  operator+ (const Factor& B) const  { return binaryOp(  B, binOpPlus()  ); };
  Factor& operator+=(const Factor& B)        { return binaryOpIP(B, binOpPlus()  ); };
  Factor  operator- (const Factor& B) const  { return binaryOp(  B, binOpMinus() ); };
  Factor& operator-=(const Factor& B)        { return binaryOpIP(B, binOpMinus() ); };
  Factor  operator* (const Factor& B) const  { return binaryOp(  B, binOpTimes() ); };
  Factor& operator*=(const Factor& B)        { return binaryOpIP(B, binOpTimes() ); };
  Factor  operator/ (const Factor& B) const  { return binaryOp(  B, binOpDivide()); };
  Factor& operator/=(const Factor& B)        { return binaryOpIP(B, binOpDivide()); };

  Factor  operator+ (const value B) const    { return binaryOp(  B, binOpPlus()  ); };
  Factor& operator+=(const value B)          { return binaryOpIP(B, binOpPlus()  ); };
  Factor  operator- (const value B) const    { return binaryOp(  B, binOpMinus() ); };
  Factor& operator-=(const value B)          { return binaryOpIP(B, binOpMinus() ); };
  Factor  operator* (const value B) const    { return binaryOp(  B, binOpTimes() ); };
  Factor& operator*=(const value B)          { return binaryOpIP(B, binOpTimes() ); };
  Factor  operator/ (const value B) const    { return binaryOp(  B, binOpDivide()); };
  Factor& operator/=(const value B)          { return binaryOpIP(B, binOpDivide()); };
  Factor  operator^ (const value B) const    { return binaryOp(  B, binOpPower() ); };
  Factor& operator^=(const value B)          { return binaryOpIP(B, binOpPower() ); };


////////////////////////////////////////////////////////////////////////////////
// Partition function, entropy, and normalization
////////////////////////////////////////////////////////////////////////////////
  Factor  normalized()  const { Factor F=*this; F.normalize(); return F;  }
  Factor& normalize() { double Z=sum(); if (Z!=0) for (vsize i=0;i<N_;i++) t_[i]/=Z; return *this; }

  //double partition()    const { double Z=0; for (vsize i=0;i<numel();i++) Z+=t_[i]; return Z; };
  double logpartition() const { return std::log(sum())/std::log(2.0); };

  value entropy(void) const { 
    value H=0, Z=0;
	for (vsize i=0;i<N_;i++) { Z += t_[i]; H -= t_[i]*std::log(t_[i]); }
	H/=Z; H+=std::log(Z);
	H/=std::log(2.0);			// base 2 log?
	return H;
  }


////////////////////////////////////////////////////////////////////////////////
// Elimination operators (sum, max, min, ...)
////////////////////////////////////////////////////////////////////////////////

  Factor sum(Variables const& sumOut) const {  
	Variables target = v_ - sumOut;
	return marginal(target);
  }
  value sum() const { value Z=0; for (vsize i=0;i<N_;i++) Z+=t_[i]; return Z; }

  Factor sumPower(Variables const& sumOut,value pow)  const {  
    if (pow==1.0)          return sum(sumOut); 
    else if(pow==-infty()) return min(sumOut);
    else if(pow== infty()) return max(sumOut); 
    else {
      Factor F=*this; F.log(); F*=pow; F=F.logsumexp(sumOut); F/=pow; F.exp();
      return F;
    }
  }

  Factor logsumexp(const Variables& sumOut) const {  
    Variables target = v_ - sumOut;
    Factor mx = maxmarginal(target);
    Factor Scaled = *this - mx;
    Scaled.exp(); 
    mx += Scaled.marginal(target).log();
    return mx;
  }

/*  // Alternative version (?)
  Factor logsumexp(const Variables& sumOut) const {  
	Variables target = v_ - sumOut;
	Factor F(target,0.0);
	subindex s(v_,target);
	for (vindex i=0;i<N_; ++i,++s) { 
	  value mx=F[s],mn=t_[i];  if (mx<mn) { mx=t_[i]; mn=F[s]; }
	  F[s]=mx+std::log(1+std::exp(mn-mx)); 
	}
	return F;
  }
*/

  Factor max(Variables const& sumOut) const {  
	Variables target = v_ - sumOut;
    return maxmarginal(target);
  }
  //value max() const { value mv=-infty(); for (vsize i=0;i<N_;i++) mv=(mv<t_[i])? t_[i] : mv; return mv; }
  value max() const { value* mx=max_element(t_,t_+N_); return *mx; }

  Factor min(Variables const& sumOut) const {  
	Variables target = v_ - sumOut;
    return minmarginal(target);
  }
  //value min() const { value mv=infty(); for (vsize i=0;i<N_;i++) mv=(mv<t_[i])? mv : t_[i]; return mv; }
  value min() const { value* mn=min_element(t_,t_+N_); return *mn; }

  inline vsize argmax() const {  
	vsize idx=0;
	for (vsize i=0;i<N_; ++i) idx=(t_[idx] > t_[i]) ? idx : i;
	return idx;
  }
  inline vsize argmin() const {  
	vsize idx=0;
	for (vsize i=0;i<N_; ++i) idx=(t_[idx] > t_[i]) ? i : idx;
	return idx;
  }

  //Factor condition(Variables const& var, std::vector<vsize> const& val) const {  }; 	//!!!

  //vsize sample() const {  };	// !!!

  Factor marginal(Variables const& target) const {  
	Factor F(target,0.0);
	subindex s(v_,target);
	for (vsize i=0;i<N_; ++i,++s) F[s]+=t_[i];
	return F;
  }
  Factor maxmarginal(Variables const& target) const {
	Factor F(target,-infty());
	subindex s(v_,target);
	for (vsize i=0;i<N_; ++i,++s) F[s]=(F[s] > t_[i]) ? F[s] : t_[i];
	return F;
  }
  Factor minmarginal(Variables const& target) const {
	Factor F(target,infty());
	subindex s(v_,target);
	for (vsize i=0;i<N_; ++i,++s) F[s]=(F[s] > t_[i]) ? t_[i] : F[s];
	return F;
  }


////////////////////////////////////////////////////////////////////////////////
// Misc other functions
////////////////////////////////////////////////////////////////////////////////
  enum DistType { DIST_L2 , DIST_L1 , DIST_LINF , DIST_KL , DIST_HPM , DIST_MAS };
  double distance(Factor const& F2, Factor::DistType type=DIST_L2) const {
    assert( vars() == F2.vars() );
    Factor F(*this);								// make a copy for manipulation
    double dist=-1.0;								// 
    value Z;
    switch (type) {
		case DIST_L2:								// L2, sum of squared errors
			F-=F2; F*=F; dist=F.sum();
			break;
		case DIST_L1: 								// L1, sum of absolute errors
			F-=F2; dist=F.abs().sum();
			break;
		case DIST_LINF:								// L-infinity, max absolute error
			F-=F2; dist=F.abs().max(); 
			break;
		case DIST_KL:								// KL-divergence (relative entropy)
			Z=sum(); F/=F2; F*=F2.sum()/Z; F.log(); F*=*this; dist=F.sum()/Z;
			break;
		case DIST_HPM:								// Hilbert's projective metric
			F/=F2; F.log(); dist=F.max()-F.min();	//   aka "dynamic range"
			break;
		case DIST_MAS:								// "MAS" error value (not a metric)
			F.log(); Factor lF2(F2); lF2.log(); F/=lF2; 
			double mx=F.max(), mn=1/F.min(); dist=((mx > mn) ? mx : mn)-1;
			break;
    }
    return dist;
  }

  double norm(Factor::DistType type=DIST_L2) const {
	Factor F(*this);								// make a copy for manipulation
	double dist=-1.0;								// 
	switch (type) {
		case DIST_L2:								// L2, sum of squared errors
			F*=F; dist=F.sum();
			break;
		case DIST_L1: 								// L1, sum of absolute errors
			dist=F.abs().sum();
			break;
		case DIST_LINF:								// L-infinity, max absolute error
			dist=F.abs().max(); 
			break;
		case DIST_KL:								// KL-divergence (relative entropy => entropy?)
			return entropy();
			break;
		case DIST_HPM:								// Hilbert's projective metric
			F.log(); dist=F.max()-F.min();			//   aka "dynamic range"
			break;
		case DIST_MAS: assert(type!=DIST_MAS); break;	// can't use this method as a norm
	}
	return dist;
  }

  // mean(F1...Fn);
  // geomean(F1...Fn);
	

  enum DecompType { DECOMP_L2 , DECOMP_L2HPM , DECOMP_L2MAS };
	std::vector<Factor> decompSum(std::vector<Variables> vlist, Factor::DecompType method) const {
	  int nF=vlist.size();
		double mx,mn;
	  std::vector<Factor> Flist(nF);

		Factor tmp,F=*this;
		switch (method) {
			case DECOMP_L2: //L2
			  double Cn,Cd;
				Cd=F.numel(); Cn=F.sum(); // /Cd*(1-1.0/nF);
				for (int j=0;j<nF;j++) {
				  Flist[j] = F.marginal( vlist[j] );
					double D = Cd/Flist[j].numel();
			 		Flist[j]/= D;
					Flist[j]-= Cn/Cd*(1.0-1.0/(nF-j));
          F  -= Flist[j];
					Cn -= Flist[j].sum()*D;
				}
				break;
			case DECOMP_L2HPM: //L2+HPM
			  Flist = decompSum(vlist,DECOMP_L2);
				for (int j=0;j<nF;j++) F-=Flist[j];
				mx=F.max(); mn=F.min();
				for (int j=0;j<nF;j++) Flist[j]+=(mx+mn)/2/nF;
				break;
			case DECOMP_L2MAS: //L2+MAS
			  Flist = decompSum(vlist,DECOMP_L2);
				F=Flist[0]; for (int j=1;j<nF;j++) F+=Flist[j];
				F /= *this; F.log(); 
				mx=F.max(); mn=F.min();
				for (int j=0;j<nF;j++) Flist[j]*=std::exp(-(mx+mn)/2/nF);
				break;
			} 
		return Flist;
	}

	std::vector<Factor> decompProd(std::vector<Variables> vlist, Factor::DecompType method) const {
		Factor F=*this; F.log();
		std::vector<Factor> Flist = F.decompSum(vlist,method);
		for (int j=0;j<vlist.size();j++) Flist[j].exp();
		return Flist;
	}
  
  void display(void) { };



  /////////////////////////////
  // Private object functions
  /////////////////////////////
 protected:
  Variables v_;				// variable list vector
  vsize  N_;
  value* t_;				// table of values
#ifdef MEX
  mxArray* mHandle;
#endif

#ifdef MEX 
  inline void mxInit() { mHandle=NULL; };
  static inline bool isfinite(value v) { return mxIsFinite(v); }
  static inline bool isnan(value v)    { return mxIsNaN(v); }
  static inline value infty()          { return mxGetInf(); }
#else
  inline void mxInit() { };
  //static inline bool isfinite(value v) { return std::abs(v)!=std::numeric_limits<value>::infinity(); }
  static inline bool isfinite(double v) { return (v <= DBL_MAX && v >= -DBL_MAX); }
  static inline bool isnan(value v)    { return (v!=v); }
  static inline value infty()          { return std::numeric_limits<value>::infinity(); }
#endif

};


////////////////////////////////////////////////////////////////////////////////
// "Static" functions that operate on Factor class variables
////////////////////////////////////////////////////////////////////////////////

inline Factor abs(const Factor& A)   { Factor F=A; F.abs(); return F; };
inline Factor exp(const Factor& A)   { Factor F=A; F.exp(); return F; };
inline Factor log(const Factor& A)   { Factor F=A; F.log(); return F; };
inline Factor log2(const Factor& A)  { Factor F=A; F.log(); F/=log(2.0); return F; };
inline Factor log10(const Factor& A) { Factor F=A; F.log10(); return F; };
//inline Factor power(const Factor& A, double pow) { Factor F=A; F.power(pow); return F; };

#endif
