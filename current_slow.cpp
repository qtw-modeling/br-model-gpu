#include "common.hpp"

// TODO: evrthg below

// gating vars
real alpha_n_CPU(real V)
{
    return 0.01 * (V + 55.0) / (1.0 - exp(-(V + 55.0) / 10.0));
}

real beta_n_CPU(real V)
{
    return 0.125 * exp(-(V + 65.0) / 80.0);
}

// S --- stands for "Slow"
real CurrentS(real V, real m, real n, real h)
{
    real ENernst = 50.;
    real gMax = 120.;

    // TODO: check sign
    return gMax * m * m * m * h * (V - ENernst);
}