#ifndef __VARIABLES_H
#define __VARIABLES_H

#include "basics.hh"
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <stdint.h>

/*
"Variables" class -- compatible with uint32 arrays of variable ids used in the matlab version
  Mainly stores variable ID list (sorted uint32) and dimensions (mwSize)
  The main C++ interface allows for set-theoretic operations (intersect, union, etc) and some calculations of interest

Some problems: creating variables without dimensions can cause problems in later functions...
Need to clean these up in the Factor class...

*/


class Variables {
 public:
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructors, Destructor 
/////////////////////////////////////////////////////////////////////////////////////////////////////////

  Variables(void); // Blank constructor
  Variables(const Variables& var); // Copy constructor from variable
  Variables(const vindex Nv, const vindex* vlist=NULL, const vsize* dlist=NULL); // Constructor from pointers
  //Variables(const vindex v); // Constructor from single index
#ifdef MEX
  Variables(mxArray* ma); // Copy constructor from Matlab vector
#endif
  // Destructor /////////////////////////
  ~Variables(void);

#ifdef MEX
  Variables& mxCopy(const mxArray* mNew, vsize* dnew);     // mxCopy : Create a *copy* of the matlab data and wrap it
  Variables& mxSet(mxArray* mNew, const vsize* dnew=NULL); // mxSet : Wrap a given matlab memory structure
  Variables& mxNew(const vsize* dnew=NULL);                // mxNew : Create new matlab memory to hold current data
  inline bool mxAvail() const { return mHandle != NULL; }; // mxAvail : Check whether we are currently a matlab matrix
  mxArray* mxGet();                                        // mxGet   : Return the matlab handle (create if !exist)
#else
  inline bool mxAvail() const { return false; };
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Accessors
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  const vindex   nvar()         const { return N_;  };   // # of variables
  const vsize*   dims()         const { return d_;  };   // dimensions of variables
  vindex& operator[] (vindex c) const { return v_[c]; }; // index into variables
  void setDimsFromGlobal(const vindex* dGlobal);         // set dims from a global vector (!!!)
  void setDimsFromGlobal(const Variables& Global);       // set dims from a global object

  bool operator== ( const Variables& B ) const; 
  bool operator!= ( const Variables& B ) const {return !(*this==B);}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copy
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  Variables& operator= (const Variables& B);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Operators:  union (+,|), setdiff (-,/), intersection (&), xor (^)
/////////////////////////////////////////////////////////////////////////////////////////////////////////
  Variables  operator+  ( const Variables& B ) const;                      // union
  Variables& operator+= ( const Variables& B ) { return (*this=*this+B); } // 
  Variables  operator|  ( const Variables& B ) const	{ return *this+B; }  // union (also)
  Variables& operator|= ( const Variables& B ) { return (*this=*this|B); } // 

  Variables  operator- (const Variables& B) const;                         // set-diff
  Variables& operator-= ( const Variables& B ) { return (*this=*this-B); } // 
  Variables  operator/ ( const Variables& B ) const	{ return *this-B; }    // set-diff (also)
  Variables& operator/= ( const Variables& B ) { return (*this=*this/B); } // 

  Variables  operator& (const Variables& B) const;                         // intersection
  Variables& operator&= ( const Variables& B ) { return (*this =(*this & B)); }	// 

  Variables  operator^ ( const Variables& B ) const;                       // set-symmetric diff (xor)
  Variables& operator^= ( const Variables& B ) { return (*this=*this^B); } // 

/* MISSING FUNCTIONS ********************

isMember
isSubset / contains / ...
add (var,dim)?
remove (var)? 
lexicographical orderings for maintaining sets of sets?
# of states? (nStates? nrStates?)

****************************************/

friend std::ostream& operator<< (std::ostream &os, Variables const& v);    // ostream output

protected:
#ifdef MEX
  mxArray* mHandle;			// matlab data structure handle (array of uint32s)
  inline bool mxTrue() { return mHandle!=NULL; }	// tests and access f'ns for matlab version
  inline void mxInit() { mHandle=NULL; }
#else
  inline bool mxTrue() { return false; }			// tests and access f'ns for non-matlab version
  inline void mxInit() { }
#endif
  vindex  N_;				// number of variables 
  vindex* v_;				// vector of variable ID #s
  const vsize*  d_;			// dimensions of the variables (potentially from Matlab)
  vsize*  dLocal;			// non-const version (equals d_) if we allocated the dimensions ourselves

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructors, Destructor 
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Blank constructor
Variables::Variables(void) { 
  N_=0; d_=dLocal=NULL; v_=NULL; mxInit();
};  

// Copy constructor from variable
Variables::Variables(const Variables& var) { 
  mxInit(); N_=var.N_; v_=NULL; d_=dLocal=NULL; if (N_) { v_=new vindex[N_]; d_=dLocal=new vsize[N_]; }
  for (vindex i=0;i<N_;i++) { v_[i]=var.v_[i]; dLocal[i]=var.d_[i]; }
}

// Constructor from pointers
Variables::Variables(const vindex Nv, const vindex* vlist, const vsize* dlist) {
  N_=Nv; mxInit(); v_=NULL; d_=dLocal=NULL;
  if (N_) { v_=new vindex[N_]; d_=dLocal=new vsize[N_]; }		// allocate with standard memory mgr
  if (vlist) for (vindex i=0;i<N_;i++) v_[i]=vlist[i];      // and copy the data off  (sort? unique? !!!)
  else       for (vindex i=0;i<N_;i++) v_[i]=0;
  if (dlist) for (vindex i=0;i<N_;i++) dLocal[i]=dlist[i];
  else       for (vindex i=0;i<N_;i++) dLocal[i]=0;
}

// Destructor /////////////////////////
Variables::~Variables(void) {
  if (!mxAvail()) {												// if it's not in MEX control,
    if (v_) delete[] v_;                 	// delete any allocated memory
  }
  if (dLocal) delete[] dLocal;    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MEX  //////////////////////////////////////////////////////////////////////////////////////////
// Copy constructor from Matlab vector
Variables::Variables(mxArray* ma) { 
  mHandle = ma; N_=mxGetN(ma); 
  v_=(vindex*) mxGetData(ma); 
  d_=dLocal=new vsize[N_];  for (vsize i=0;i<N_;i++) dLocal[i]=0;
};

// mxCopy : Create a *copy* of the matlab data and wrap it
Variables::Variables& mxCopy(const mxArray* mNew, vsize* dnew) {
  //std::cout<<"Copying variable for wrap\n";
  mxArray* mCopy = mxDuplicateArray(mNew);                            // duplicate the data (deep copy)
  mxSet(mCopy);                                                       //   and wrap the resulting matrix
};

// mxSet : Wrap a given matlab memory structure
Variables::Variables& mxSet(mxArray* mNew, const vsize* dnew=NULL) { 
  //std::cout<<"Setting variable wrapper, "<<mxAvail()<<"\n";
  if (mxAvail()) mxDestroyArray(mHandle);                             // destroy any old data we were storing
  else if (v_) delete[] v_;                                           //  ""
  mHandle=mNew;                                                       // save our new structure
  N_=mxGetN(mNew);                                                    // get the new data we will store, size
  v_ = (vindex*) mxGetData(mNew);                                     //   and values
  if (dnew) { d_=dnew; if (dLocal) { delete[] dLocal; dLocal=NULL; }}	// take new dims if they're passed
  else      { d_=dLocal=new vsize[N_]; for (vsize i=0;i<N_;i++) dLocal[i]=0; }  // otherwise create "blanks"
  if (N_==0) {v_=NULL; d_=dLocal=NULL; }
  return *this;
};

// mxNew : Create a new matlab memory structure to hold our current data
Variables::Variables& mxNew(const vsize* dnew=NULL) {
  //std::cout<<"Creating new variable wrapper\n";
  mxArray* mNew = mxCreateNumericMatrix(1,N_,mxUINT32_CLASS,mxREAL);   // create new memory (a matrix)
  vindex*  vNew = (vindex*) mxGetData(mNew);                           //  and get its data pointer
  for (vindex i=0;i<N_;i++) vNew[i]=v_[i];                             // copy off our current data
  if (dnew) { d_=dnew; if (dLocal) { delete[] dLocal; dLocal=NULL; }}  // take new dims if they're passed
  if (mxAvail()) {                                                     // if we were already a matlab matrix
    mxDestroyArray(mHandle);                                           //  destroy our old version
  } else {                                                             // if we weren't, de-allocate the
    if (v_) delete[] v_;                                               //  memory we were using
  }
  v_=vNew; mHandle=mNew;                                               // save our new structure
  return *this;
};

// mxGet   : Return current matlab handle (create if doesn't exist)
Variables::mxArray* mxGet() {
  if (!mxAvail()) mxNew();
  return mHandle;
};
#endif  //////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Accessors
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Variables::setDimsFromGlobal(const Variables& Global) { 
  assert(dLocal!=NULL); 
  Variables tmp = Global & *this;
  for (vindex i=0;i<N_;++i) dLocal[i]=tmp.d_[i];
};
void Variables::setDimsFromGlobal(const vindex* Global) { 
  assert(dLocal!=NULL); 
  for (vindex i=0;i<N_;++i) dLocal[i]=Global[v_[i]];
};

bool Variables::operator== ( const Variables& B ) const {
  if (N_ != B.N_) return false;
  for (vindex i=0;i<N_;++i) if (v_[i]!=B.v_[i]) return false;
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copy
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Variables& Variables::operator= (const Variables& B) {
  if (N_>=B.N_) { N_=B.N_; for (vindex i=0;i<N_;++i) v_[i]=B.v_[i]; }		// if we can store it as is, do so
  else {                                                                // otherwise, delete, allocate, & copy
    N_=B.N_;
    assert(!mxAvail());                                                 // better not be a matlab object...
    if (v_) { delete[] v_; }          v_=new vindex[N_];
    if (dLocal) { delete[] dLocal; } 	d_=dLocal=new vsize[N_];
    if (N_==0) {v_=NULL; d_=dLocal=NULL; }
    for (vindex i=0;i<N_;++i) { v_[i]=B.v_[i]; dLocal[i]=B.d_[i]; }
  }
  return *this;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Operators:  union (+,|), setdiff (-,/), intersection (&), xor (^)
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Variables Variables::operator+ ( const Variables& B ) const {
  Variables dest(N_+B.N_);
  vindex i,j,k;
  for (i=0,j=0,k=0;i<N_&&j<B.N_;) {
    if (v_[i]==B[j])     { dest.dLocal[k]=  d_[i]; dest.v_[k++]=  v_[i++]; j++; }
	  else if (v_[i]>B[j]) { dest.dLocal[k]=B.d_[j]; dest.v_[k++]=B.v_[j++];      }
    else                 { dest.dLocal[k]=  d_[i]; dest.v_[k++]=  v_[i++];      }
  }
  while (i<N_)   { dest.dLocal[k]=  d_[i]; dest[k++]=  v_[i++]; }
  while (j<B.N_) { dest.dLocal[k]=B.d_[j]; dest[k++]=B.v_[j++]; }
  dest.N_=k;
  return dest;
};

Variables Variables::operator- (const Variables& B) const {
  Variables dest(N_);
  vindex i,j,k;
  for (i=0,j=0,k=0;i<N_&&j<B.N_;) {
    if      (v_[i]< B[j])  { dest.dLocal[k]=d_[i]; dest.v_[k++]=v_[i++]; }
    else if (v_[i]==B[j])  { i++; j++; }
    else                   { j++;      }
  }
  while (i<N_) { dest.dLocal[k]=d_[i]; dest.v_[k++]=v_[i++]; }
  dest.N_=k;
  return dest;
}

Variables  Variables::operator& (const Variables& B) const {
  Variables dest( N_>B.N_ ? B.N_ : N_ );
  vindex i,j,k;
  for (i=0,j=0,k=0;i<N_&&j<B.N_;) {
    if      (v_[i]< B[j]) { i++; }
    else if (v_[i]==B[j]) { dest.dLocal[k]=d_[i]; dest.v_[k++]=v_[i++]; j++; }
    else                  { j++; }
  }
	dest.N_=k;
  return dest;
}

Variables  Variables::operator^ ( const Variables& B ) const {
  Variables dest(N_+B.N_);
  vindex i,j,k;
  for (i=0,j=0,k=0;i<N_&&j<B.N_;) {
    if (v_[i]==B[j])     { i++; j++; }
    else if (v_[i]>B[j]) { dest.dLocal[k]=B.d_[j]; dest.v_[k++] = B.v_[j++];      }
    else                 { dest.dLocal[k]=  d_[i]; dest.v_[k++] =   v_[i++];      }
  }
  while (i <   N_)        { dest.dLocal[k]=  d_[i]; dest.v_[k++] =   v_[i++]; }
  while (j < B.N_)        { dest.dLocal[k]=B.d_[j]; dest.v_[k++] = B.v_[j++];  }
  dest.N_ = k;
  return dest;
};


std::ostream& operator<< (std::ostream &os, Variables const& v) {
  os<<"{";
  for (vindex i=0;i<v.nvar();i++) os << ( (i) ? "," : "" ) << v[i];
  os<<"}";
	return os;
}

#endif 
