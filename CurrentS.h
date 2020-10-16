#ifndef __CURRENT_S_H
#define __CURRENT_S_H

#include "common.h"

#pragma acc routine
real alpha_d(real V);
#pragma acc routine
real beta_d(real V);

#pragma acc routine
real alpha_f(real V);
#pragma acc routine
real beta_f(real V);

#pragma acc routine
real d_inf(real V);
#pragma acc routine
real f_inf(real V);

// MAIN FUNC: CURRENT
#pragma acc routine
real CurrentS(real V, real d, real f, real concCa);

#endif