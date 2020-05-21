#include "common.hpp"

// GPU vers of funcs
#pragma acc routine
real alpha_x(real V)
{
    //return 9. * 1e-3 / (1. + exp( -(V + 3.8)/9.71)) + 6e-4;

    return AlphaGeneralForm(V, 0.0005, 0.083, 50., 0., 0., 0.057, 1.);
}


#pragma acc routine
real beta_x(real V)
{
    //return 2.25 * 1e-4 * (V + 40.) / ( exp((V + 40.) / 13.3) - 1 ) ;

    return AlphaGeneralForm(V, 0.0013, -0.06, 20., 0., 0., -0.04, 1.);
}



#pragma acc routine
real x_inf(real V)
{
    return alpha_x(V) / (alpha_x(V) + beta_x(V));
}


#pragma acc routine
real CurrentX(real V, real x)
{
    // another K-current
    
    return 0.8*x* ( exp(0.04*(V + 77.)) - 1. ) / ( exp(0.04*(V + 35.)) );
}