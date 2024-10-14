
#ifndef HEADER_P
#define HEADER_P

#include <math.h>

// parameters for homomorphic tallying
const long N           = 4096;
const long q           = 4611686018326724609;
const long p           = 10000009;
const long tau         = 1;
const long lambda      = 128;
const long B           = 1;
const double MLinInv   = 1/sqrt(3);
const double MAInv     = 1/3.0;
const double Sigma_A   = 636904594;
const double Bound_A   = 73;    // # bits * 2
const double TwoSigmaA2 = 2*Sigma_A*Sigma_A;
const double Sigma_C   = 43238;
const double Bound_C   = 44.8;  // # bits * 2
const double TwoSigmaC2 = 2*Sigma_C*Sigma_C;

// parameters for mix-net based voting
// const long N           = 4096;
// const double q         = 4611686018326724609;
// const long p           = 2;
// const long tau         = 4096;
// const long lambda      = 128;
// const long B           = 1;
// const double MLinInv   = 1/sqrt(3);
// const double MAInv     = 1/3.0;
// const double Sigma_A   = 90675893177;
// const double Bound_A   = 88;    // # bits * 2
// const double TwoSigmaA2 = 2*Sigma_A*Sigma_A;
// const double Sigma_C   = 43238;
// const double Bound_C   = 44.8;  // # bits * 2
// const double TwoSigmaC2 = 2*Sigma_C*Sigma_C;


#endif
