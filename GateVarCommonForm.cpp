#include "common.hpp"

real AlphaGeneralForm(real Vm,
                      real c1, real c2, real c3, real c4, real c5,
                      real c6, real c7)
{
    // Function, containing a common formula for gating vars' rate funcs "alpha_i, beta_i";
    // contains 5 undetemined constants-parameters {(c0), c1, c2, c3, c4, c5} and {V0};
    // we begin with i = 1, ..., 5 (0 is excluded for correspondance with article abbreviations)

    return (c1*exp(c2*(Vm + c3)) + c4*(Vm + c5)) / (exp(c6*(Vm + c3)) + c7);
}