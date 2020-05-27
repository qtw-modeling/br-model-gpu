#pragma once

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
//#include <fstream>
//#include <cstdlib>
//#include <iostream>

typedef double real;

#define M_PI acos(-1.)

//#define DEBUG() printf("DEBUG pause\n");


// gating vars func
#pragma acc routine
real AlphaGeneralForm(real Vm,
                      real c1, real c2, real c3, real c4, real c5,
                      real c6, real c7);

// Currents are in the end of enum: "concCa" is a state var, and should be before Currents for correct enum numbering usage
//enum vars {V_, m_, h_, J_, d_, f_, x_, concCa_, INa_, IK_, IX_, IS_};