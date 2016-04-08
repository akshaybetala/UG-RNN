// Mode is: -*- C++ -*- 
// --- General.h --- Created: Sun Sep 21 16:10:50 1997
// 
// Time-stamp: <October 27, 1997 17:13:27 paolo@tequila.dsi.unifi.it>
// Update Count: 10
// 
// RCS: $Id$   $Locker$   $Log$

// = FILENAME
//    General.h
// = LIBRARY
// 
// = AUTHOR
// Paolo Frasconi (paolo@mcculloch.ing.unifi.it)
// = COPYRIGHT
// Copyright (C) 1997 Paolo Frasconi

#ifndef General_h
#define General_h 1
#define Float double

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <time.h>
//#include <values.h>

using namespace std;
#ifdef WANT_ERRNO
#include <errno.h>
#define FAULT(s) {cerr << "Error: in function " << __PRETTY_FUNCTION__ \
  << ", file " << __FILE__ << " line " << __LINE__ << ": " << s << "\n" \
  << "errno= " << errno << " (" << strerror(errno) << ")\n"; abort();}

#else
#define FAULT(s) {cerr << "Error: in function " << __PRETTY_FUNCTION__ \
  << ", file " << __FILE__ << " line " << __LINE__ << ": " << s << "\n" \
  ; abort();}
#endif
#include <string.h>
// Customization

#define CPP_COMMAND "./cpp -P -undef "

#define GZCAT_COMMAND "./zcat"

// end customization

extern int debugLevel;
#ifndef __GNUG__
#define __PRETTY_FUNCTION__ "unknown()"
#endif

#define WARN(level, s) {if (debugLevel>=level) {cerr << "Warning: in function " << __PRETTY_FUNCTION__ << ", file " << __FILE__ << " line " << __LINE__ << ": " << s << "\n";}}

#define MESSAGE(s) {cerr << "Notice: in " << __PRETTY_FUNCTION__ << ": " << s << "\n";}



#define DEBUG(level, s) {if (debugLevel >= level) {cerr << s << "\n";}}

#ifdef CHECKING
#define CHECK(condition, errorMsg) if ((condition)) {FAULT(errorMsg);}
#else
#define CHECK(condition, errorMsg) ;
#endif


#define UNKNOWN_REAL_VAL MINFLOAT
#define UNKNOWN_NOMINAL_VAL MINFLOAT
#define UNKNOWN_VAL MINFLOAT

int argmax(const Float* x, int n);
void meanvar(Float* v, int n, Float& mean, Float& std);
void printvec(const Float* x, int n);
double get_CPU_time_usage(clock_t clock1,clock_t clock2);


#endif // General_h
