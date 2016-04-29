#include <iostream>

class subindex {
  protected:
  public:
   vsize  idx_,end_;
   vsize  Nd;
   vsize  *state;
   const vsize *dims;
   bool *skipped;
   vsize *add;
   vsize *subtract;
   

   subindex(const Variables& full, const Variables& sub) {
	 	idx_=0; 
		end_=1;
		Nd=full.nvar(); 
	 	dims     = full.dims(); 
	 	state    = new vsize[Nd];
	 	add      = new vsize[Nd];
	 	subtract = new vsize[Nd];
	 	skipped  = new bool[Nd];
	 	const vsize *subd=sub.dims(), *fulld=full.dims();
	 	// Compute reference index updates
	 	vsize i,j; 
	 	for (i=0,j=0;i<Nd;++i) {
	 	 state[i]=1; 
	   skipped[i] = ( j>=sub.nvar() || sub[j]!=full[i] );
	   if (i==0) add[i]=1; else add[i]=add[i-1]*(skipped[i-1]?1:dims[i-1]);
	   subtract[i]=add[i]*((skipped[i]?1:dims[i])-1);
	   if (!skipped[i]) j++;
		 end_ *= dims[i];
	 	} 
		end_--;
   }

   ~subindex(void) {
	 delete[] state; delete[] skipped; delete[] add; delete[] subtract;
   }

	 subindex& reset() {
		 for (vindex i=0;i<Nd;++i) state[i]=1;
		 idx_=0;
		 return *this;
	 };

	 size_t end() { return end_; }

   subindex& operator++ (void) {											// prefix addition
	 for (vindex i=0;i<Nd;++i) {
	   if (state[i]==dims[i]) {
		 	state[i]=1;
		 	if (!skipped[i]) idx_ -= subtract[i];
	   } else {
		 	++state[i];
		 	if (!skipped[i]) idx_ += add[i];
		 	break;
	   }
	 }
	 return *this;
  }
  void operator++ (int) { ++(*this); }										// postfix addition (no return !!!)

  operator size_t() const { return idx_; };

};

