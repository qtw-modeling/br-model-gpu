#ifndef __CURRENT_X_H
#define __CURRENT_X_H

#include "common.h"

#pragma acc routine
real alpha_x(real V);
#pragma acc routine
real beta_x(real V);

#pragma acc routine
real x_inf(real V);

#pragma acc routine
real CurrentX(real V, real x);

#endif