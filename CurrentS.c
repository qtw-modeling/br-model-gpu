#include "common.h"

#pragma acc routine
real alpha_d(real V)
{
    // new gate var (if compared to YNI-model)

    return AlphaGeneralForm(V, 0.095, -0.01, -5., 0., 0., -0.072, 1.);
}

#pragma acc routine
real beta_d(real V)
{
    // new gate var (if compared to YNI-model)
    return AlphaGeneralForm(V, 0.07, -0.017, 44., 0., 0., 0.05, 1.);
}

#pragma acc routine
real alpha_f(real V)
{
    return AlphaGeneralForm(V, 0.012, -0.008, 28., 0., 0., 0.15, 1.);
}

#pragma acc routine
real beta_f(real V)
{
    return AlphaGeneralForm(V, 0.0065, -0.02, 30., 0., 0., -0.2, 1.);
}

#pragma acc routine
real d_inf(real V)
{
    return alpha_d(V) / (alpha_d(V) + beta_d(V));
}

#pragma acc routine
real f_inf(real V)
{
    return alpha_f(V) / (alpha_f(V) + beta_f(V));
}

// MAIN FUNC: CURRENT
#pragma acc routine
real CurrentS(real V, real d, real f, real concCa)
{
    // Calcium current, IS

    return 0.09 * d * f * (V + 82.3  + 13.0287*log(concCa) );
}