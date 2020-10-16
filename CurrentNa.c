#include "common.h"

// GPU vers of funcs
// commented return --- means dat the func has been modded from YNI to BR model


#pragma acc routine
real alpha_m(real V)
{
    //return (V + 37.) / (1. - exp( -(V + 37.)/10. )); 

    // passing coeffs from the table for BL-model in "Computing electr. act. of the heart"
    //return AlphaGeneralForm(V, 0., 0., 47., -1., 47., -0.1, -1.);

    // BR-DR model formulation
    return 0.9*(V + 42.65) / (1 - exp(-0.22*(V + 42.65))); // is it "+" or "-" in the article? TODO: check
}


#pragma acc routine
real beta_m(real V)
{
    //return 40. * exp(-5.6 * 1e-2 * (V + 62.));

    //return AlphaGeneralForm(V, 40., -0.056, 72., 0., 0., 0., 0.);

    // BR-DR model formulation
    return 1.437 * exp(-0.085*(V + 39.75)); // is it "5" or "8" in the article? TODO: check

}

#pragma acc routine
real alpha_h(real V)
{
    //return 1.209 * 1e-3 * exp( -(V + 20.)/6.534);

    //return AlphaGeneralForm(V, 0.126, -0.25, 77., 0., 0., 0., 0.);

    // BR-DR model formulation
    return 0.1 * exp( -0.193*(V + 79.65) );

}

#pragma acc routine
real beta_h(real V)
{
    //return 1. / ( exp(-(V + 30.)/10.) + 1.);

    //return AlphaGeneralForm(V, 1.7, 0., 22.5, 0., 0., -0.082, 1.);

    // BR-DR model formulation
    return 1.7 / (1 + exp(-0.095*(V + 20.5)));
}

#pragma acc routine
real alpha_j(real V)
{
    // new gate var (if compared to YNI-model)

    return AlphaGeneralForm(V, 0.055, -0.25, 78., 0., 0., -0.2, 1.);

}

#pragma acc routine
real beta_j(real V)
{
    // new gate var (if compared to YNI-model)

    return AlphaGeneralForm(V, 0.3, 0., 32., 0., 0., -0.1, 1.);
}


#pragma acc routine
real m_inf(real V)
{
    return alpha_m(V) / (alpha_m(V) + beta_m(V));
}


#pragma acc routine
real h_inf(real V)
{
    return alpha_h(V) / (alpha_h(V) + beta_h(V));
}

#pragma acc routine
real j_inf(real V)
{
    // new gate var (if compared to YNI-model)
    return alpha_j(V) / (alpha_j(V) + beta_j(V));
}


#pragma acc routine
real CurrentNa(real V, real m, real h, real J)
{
    //real ENernst = 30.;
    //real gMax = 120.;

    // BR-DR model formulation (Drouhrard-Roberge, 1987);
    // in this formulation, there is no "J"-gating variable
    return 15.*m*m*m*h*(V - 40.);

    //return (4.* m * m * m * h * J + /* NaCa exchanger*/ 0.003) * (V - 50.); // formula from original BR-model article

    // gNa = 23 (the val. described in A.Tveito's article on FDM for BL-model ) instead of old val. gNa = 4;
    //return (/* 23. */ 10. * m * m * m * h * J + /* NaCa exchanger*/ 0.003) * (V - 50.);
}