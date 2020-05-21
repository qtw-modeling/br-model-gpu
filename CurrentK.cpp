#include "common.hpp"

real CurrentK(real V)
{

    return 0.35* ( 
        4.*( exp(0.04*(V + 85.)) - 1. ) / ( exp(0.08*(V + 53.)) + exp(0.04*(V + 53.)) )
        + 0.2*(V + 23.) / (1. - exp(-0.04*(V + 23.))) // it is maybe, wrong "tail" in the formula
        );
}