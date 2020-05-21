#pragma once

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <cstdlib>
#include <iostream>

typedef double real;

#define DEBUG() printf("DEBUG pause\n"); std::cin.get()


// gating vars func
real AlphaGeneralForm(real Vm,
                      real c1, real c2, real c3, real c4, real c5,
                      real c6, real c7);