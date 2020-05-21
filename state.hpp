#include "common.hpp"

// main (independent) variables 
struct State {
    real V, m, h, J, d, f, x;
};

// extra (dependent) variables
// TODO: add at later modification of the prorgam

struct Currents {
    real* INa, *IK, *IX, *IS;
};

struct Concentrations {
    real Ca;
};