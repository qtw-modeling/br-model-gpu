//
// Created by QuoTheWhite on 27/03/2019.
//
//#include <stdio.h>
//#include <math.h>
//#include <time.h>
//#include <fstream>
//#include <cstdlib>
//#include "openacc.h"


//#include <string>
//#include <map>
//#include <cmath>
//#include <math.h> // 4some specific constants

// 4 C++11 includes
//#include <vector>

//typedef double real;
#include "common.h" // one for all includes
#include "state.h"
#include "CurrentNa.h"
#include "CurrentK.h"
#include "CurrentX.h"
#include "CurrentS.h"

// grid parameters
//#define SINGLE_CELL_MODEL

//#define numSegmentsX 20 //10 // 31x31 cell --- is not enoght for reentry to "live" in the tussie sample; more cells are needed
//#define numSegmentsY 20 //10
//#define numPointsX (numSegmentsX + 1) // = numCells
//#define numPointsY (numSegmentsY + 1) // = numCells
//#define numPointsTotal (numPointsX * numPointsY)
#define hx 0.07 //1. // uncomment if cells are connected // (1./numSegmentsX)
#define hy 0.07 //1. // uncomment if cells are connected // (1./numSegmentsY)
#define T (30.) //(1000.) // old val: 500 // endtime (in ms)
#define dt 0.005 // old val = 1e-4 // timestep (in ms)

// model parameters
#define Cm 1.
#define VRest (-80.) // NOTE: there exists no resting potential for SA node

#define APD0 (330.) // in (ms) // from Aliev-Panfilov model
#define LongitudeStim (250.) /* the value (520) results in "smoothest" action potential */ // in (ms)
#define PeriodStim (LongitudeStim * 2.) // in (ms) //APD0
#define c0 (0.01*1e-7) // concentration initial

// tissue parameters
#define Dx (50e-3) // 30e-3 --- ok for reentry
#define Dy (50e-3) // 30e-3 --- ok for reentry

// Currents are in the end of enum: "concCa" is a state var, and should be before Currents for correct enum numbering usage
enum vars {V_, m_, h_, J_, d_, f_, x_, concCa_, INa_, IK_, IX_, IS_};




real* MemAlloc(int n)
{
    return (real*)malloc(n * sizeof(real));
}


void Write2VTK(const int n, real* p, const real h, const int step)
{
    char fn[256];
    sprintf(fn, "./output/V.%d.vtk", step); 

    // C code
    FILE* f = fopen(fn, "w");
    fprintf(f, "# vtk DataFile Version 3.0\nSolution\nASCII\nDATASET RECTILINEAR_GRID\n");
    fprintf(f, "DIMENSIONS %d %d 1\n", n + 1, n + 1);
    fprintf(f, "X_COORDINATES %d double\n", n + 1);

    int i;
    for (i = 0; i < n +1; i++)
    {
        fprintf(f, "%.2e ", i*h);
    }
    fprintf(f, "\n");

    fprintf(f, "Y_COORDINATES %d double\n", n + 1);
    for (i = 0; i < n +1; i++)
    {
        fprintf(f, "%.2e ", i*h);
    }
    fprintf(f, "\n");

    fprintf(f, "Z_COORDINATES 1 double\n0\n");
    fprintf(f, "CELL_DATA %d\n", (n * n));
    fprintf(f, "SCALARS V_membrane double\nLOOKUP_TABLE default\n");    
    
    int j;
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++)
            fprintf(f, "%.2e ",  p[j * n + i]);
        fprintf(f, "\n");
    }

    fclose(f);


/* C++ code
    std::fstream f(fn, std::ios::out);
    f << "# vtk DataFile Version 3.0" << std::endl;
    f << "Solution" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET RECTILINEAR_GRID" << std::endl;
    f << "DIMENSIONS " << n + 1 << " " << n + 1 << " 1" << std::endl;
    f << "X_COORDINATES " << n + 1 << " double" << std::endl;
    for (int i = 0; i < n + 1; i++)
        f << i * h << " ";
    f << std::endl;
    f << "Y_COORDINATES " << n + 1 << " double" << std::endl;
    for (int i = 0; i < n + 1; i++)
        f << i * h << " ";
    f << std::endl;
    f << "Z_COORDINATES 1 double\n0" << std::endl;
    f << "CELL_DATA " << (n * n) << std::endl;
    f << "SCALARS V_membrane double\nLOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++)
            f << p[j * n + i] << " ";
        f << std::endl;
    }
    f.close();
    */
}


int CalculateLinearCoordinate_CPU(int i, int j, int numPointsX) {
    return i + j*numPointsX;
}

#pragma acc routine
int CalculateLinearCoordinate(int i, int j, int numPointsX) {
    return i + j*numPointsX;
}



real FuncSinglePeriod(real t)
{
    // We consider onlt t > 0! //

    if (t < PeriodStim / 2.)
        return 0.; // first: no stimulation
    else
        return 1.; // second: stimulation


/*
 // "__П__ " --- func looks like dis
    if (t < 2./7.*PeriodStim || t > 4./7*PeriodStim)
        return 0.;
    else 
        return 1.;
*/

}


real IsStim(real t)
{

    // We consider only t > 0! //

    // Periodic stimulation, with each stimulation's duration = LongitudeStim
    real offset = (int)(t / PeriodStim) * PeriodStim;
    return FuncSinglePeriod(t - offset);


/*
    // permanent stimulation
    if (t > 400) // times in (ms): let 100 sec for transient, after this --- stimulate once
        return 1.; // stim
    else
        return 0.; // no stim
*/
    



    //return 1.; // stimulus is always ON
}


#pragma acc routine
real I_Stim(int i, int j, real value)
{
    //int x0 = (int)((real)numSegmentsX /2.); int y0 = (int)((real)numSegmentsY /2.); // tmp values

    /* uncomment if cells are connected
    if ((i > 1 && i < numPointsX - 2) && (j == 2))
        return value;
    else
        return 0.;
    */

    return value;
}

real TotalIonCurrent_CPU(int idx, real V, real m, real h, real j, real d, real f, real x,
                     real *INa, real *IK, real *IX, real *IS, real concCa)
{
    INa[idx] = CurrentNa(V, m, h, j);
    IK[idx] = CurrentK(V);
    IX[idx] = CurrentX(V, x);
    IS[idx] = CurrentS(V, d, f, concCa);

    // TODO: check the sign of the expression below
    return -(INa[idx] + IK[idx] + IX[idx] + IS[idx]);
}


#pragma acc routine
real TotalIonCurrent(int idx, real V, real m, real h, real j, real d, real f, real x,
                     real *INa, real *IK, real *IX, real *IS, real concCa)
{
    INa[idx] = CurrentNa(V, m, h, j);
    IK[idx] = CurrentK(V);
    IX[idx] = CurrentX(V, x);
    IS[idx] = CurrentS(V, d, f, concCa);

    // TODO: check the sign of the expression below
    return -(INa[idx] + IK[idx] + IX[idx] + IS[idx]);
}




real TimeViaPhase(real phase_, real Period_, real t0_) 
{
    /* phi = phase; Period = T; t0 = time of first depolarization. */
    // underscores --- for not confusing with global scope area names (or #defines)
    return t0_ + Period_ / (2.*M_PI) * phase_;
    
}

// phase calculations
real* VviaPhase(real phase) 
{
    // we dont need to perform memalloc for members within these structs
    struct State* Old = (struct State*)malloc(sizeof(struct State));
    struct State* New = (struct State*)malloc(sizeof(struct State));
    struct Concentrations* ConcentrationsOld = (struct Concentrations*)malloc(sizeof(struct Concentrations));
    struct Concentrations* ConcentrationsNew = (struct Concentrations*)malloc(sizeof(struct Concentrations)); //new Concentrations;

    // initial conditions: only for Old structs: New will be calculated in the loop
    Old->V = VRest;
    Old->m = 0.5;                //m_inf_CPU(VRest); // 0.5
    Old->h = 0.5; //h_inf_CPU(VRest); // 0.5
    Old->J = 0.5;    // 0.1; // may be false; TODO perform calcs with higher T
    Old->d = 0.5; //0.;
    Old->f = 0.5; //1.;
    Old->x = 0.5; // TODO: check this "initial condition" to result in a limit cycle
    ConcentrationsOld->Ca = c0; // 0.1*1e-7; // random value

    //std::map<std::string, real> stateOfPhase;
    real* stateOfPhase = MemAlloc(8); // 8 --- number of vars., wiht currents excluded
#ifdef SINGLE_CELL_MODEL
    stateOfPhase["V"] = Old->V;
    stateOfPhase["m"] = Old->m;
    stateOfPhase["h"] = Old->h;
    stateOfPhase["J"] = Old->J;
    stateOfPhase["d"] = Old->d;
    stateOfPhase["f"] = Old->f;
    stateOfPhase["x"] = Old->x;
    stateOfPhase["concCa"] = ConcentrationsOld->Ca;
    return stateOfPhase;
#endif


    // ...but here --- YES
    struct Currents* CurrentsOld = (struct Currents*)malloc(sizeof(struct Currents)); //new Currents;
    
    // ??? maybe we dont need CurrentsNew??? we are not swapping OLD and NEW
    struct Currents* CurrentsNew = (struct Currents*)malloc(sizeof(struct Currents));
    // dont forget to memalloc pointers within Currents struct;
    // we use only single cell
    CurrentsOld->INa = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsOld->IK = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsOld->IX = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsOld->IS = (real*)malloc(sizeof(real)); //new real[1];
    
    // ...same for "New" struct here
    CurrentsNew->INa = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsNew->IK = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsNew->IX = (real*)malloc(sizeof(real)); //new real[1];
    CurrentsNew->IS = (real*)malloc(sizeof(real)); //new real[1];



    real THRESHOLD = -30.; // hardcoded for now
    real tEndTransient = 100.; // in (ms)
    real time0; // random value
    real time1; // random value
    bool isThresholdFound = false;

    real tCurrent = 0;
    int counter = 0;
    // TODO:
    // main loop: timestepping
    while (1)
    {

        ///////////////// gating variables: ode ("reaction") step
        // TODO: make only ONE read of Old->V, etc. from memory; to more speedup, esp. for GPU
        New->m = m_inf(Old->V) + (Old->m - m_inf(Old->V)) * exp(-dt * (alpha_m(Old->V) + beta_m(Old->V)));

                    
        New->h = h_inf(Old->V) + (Old->h - h_inf(Old->V)) * exp(-dt * (alpha_h(Old->V) + beta_h(Old->V)));
        
        New->J = j_inf(Old->V) + (Old->J - j_inf(Old->V)) * exp(-dt * (alpha_j(Old->V) + beta_j(Old->V)));

        New->d = d_inf(Old->V) + (Old->d - d_inf(Old->V)) * exp(-dt * (alpha_d(Old->V) + beta_d(Old->V)));

        New->f = f_inf(Old->V) + (Old->f - f_inf(Old->V)) * exp(-dt * (alpha_f(Old->V) + beta_f(Old->V)));

        New->x = x_inf(Old->V) + (Old->x - x_inf(Old->V)) * exp(-dt * (alpha_x(Old->V) + beta_x(Old->V)));

        // membrane potential calc                  
        New->V = Old->V + dt / Cm * ( TotalIonCurrent_CPU(0 /* 0 --- index, stands for the single cell model */ , 
                                 Old->V, Old->m, Old->h, Old->J, Old->d, Old->f, Old->x,
                                 (CurrentsOld->INa), (CurrentsOld->IK), (CurrentsOld->IX), 
                                 (CurrentsOld->IS), ConcentrationsOld->Ca)
                                 + IsStim(tCurrent)*I_Stim(0, 0, 2.1) ); // "standart" I_stim = 1e

        // concentrations calc
        ConcentrationsNew->Ca = ConcentrationsOld->Ca + dt * (
                                                            -1e-7*CurrentsOld->IS[0] 
                                                            + 0.07*(1e-7 - ConcentrationsOld->Ca) 
                                                            ); // index "0" --- for array of length 1
        
        // when threshold is found
        if ((tCurrent >= tEndTransient) && (New->V > THRESHOLD) && (Old->V < THRESHOLD))
        {
            // when 2nd threshold time (t1) is found: set t1 and then exit the loop
            if (isThresholdFound == true)
            {

                //printf("Iteration #%d; V_membrane: %.2f\n", counter, Old->V);
                //DEBUG();
                time1 = tCurrent; // nearest-neighbour interpolaion; change to linear!
                break;            // phase(V)
            }
            else // when threshold time (t0) is found: set t0
            {
                //printf("Iteration #%d; V_membrane: %.2f\n", counter, Old->V);
                //DEBUG();
                time0 = tCurrent; // nearest-neighbour interpolaion; change to linear!
                isThresholdFound = true;
                //return ; // phase(V)
            }

        }

        
        tCurrent += dt;
        counter += 1;

        // swapping time layers
        struct State* tmp;
        tmp = Old; Old = New; New = tmp;

        struct Concentrations* tmpConc; // for swapping CONCENTRATIONS-type structs
        tmpConc = ConcentrationsOld; ConcentrationsOld = ConcentrationsNew; ConcentrationsNew = tmpConc;

        //printf("Iteration #%d; V_membrane: %.2f\n", counter, Old->V);
        //DEBUG();

    } // while

    //printf("t0 = %.2f, t1 = %.2f\n", time0, time1);
    //std::cin.get();

    // set vars, calculated within the loop
    real period = (time1 - time0);//*0.5; // period of oscillations; remove "0.5" when period calc bug is found!
    //printf("First loop is finished; period of oscillations: %.2f ms\n", period);
    
    
    // repeat the loop (calculations) again and find V(phi)
    tCurrent = 0; // again
    
    real tOfPhase = TimeViaPhase(phase, period, time0);
    //printf("Phase: %.2f, tOfPhase: %.2f\n", phase, tOfPhase);
    //std::cin.get();

    //real VOfPhase; // to be determined in the loop below
    //std::map<std::string, real> stateOfPhase;

    // (again) initial conditions: only for Old structs: New will be calculated in the loop
    Old->V = VRest;
    Old->m = 0.5; //m_inf_CPU(VRest); // 0.5
    Old->h = 0.5; //h_inf_CPU(VRest); // 0.5
    Old->J = 0.5;    // 0.1; // may be false; TODO perform calcs with higher T
    Old->d = 0.5;    //0.;
    Old->f = 0.5;    //1.;
    Old->x = 0.5;    // TODO: check this "initial condition" to result in a limit cycle

    ConcentrationsOld->Ca = c0; // rand value


    // (again): main loop: timestepping
    while (1)
    {
        // it means, dat we found the moment of time, corresponding to the phase value
        if (tCurrent >= tOfPhase)
        {    
            //VOfPhase = Old->V; // nearest-neighbour iterpolation; change to linear!
            stateOfPhase[V_] = Old->V;
            stateOfPhase[m_] = Old->m;
            stateOfPhase[h_] = Old->h;
            stateOfPhase[J_] = Old->J;
            stateOfPhase[d_] = Old->d;
            stateOfPhase[f_] = Old->f;
            stateOfPhase[x_] = Old->x;
            stateOfPhase[concCa_] = ConcentrationsOld->Ca;
            
            break; // exit the loop when found required t of phase
            
        }

        ///////////////// gating variables: ode ("reaction") step
        // TODO: make only ONE read of Old->V, etc. from memory; to more speedup, esp. for GPU
        New->m = m_inf(Old->V) + (Old->m - m_inf(Old->V)) * exp(-dt * (alpha_m(Old->V) + beta_m(Old->V)));

                    
        New->h = h_inf(Old->V) + (Old->h - h_inf(Old->V)) * exp(-dt * (alpha_h(Old->V) + beta_h(Old->V)));
        
        New->J = j_inf(Old->V) + (Old->J - j_inf(Old->V)) * exp(-dt * (alpha_j(Old->V) + beta_j(Old->V)));

        New->d = d_inf(Old->V) + (Old->d - d_inf(Old->V)) * exp(-dt * (alpha_d(Old->V) + beta_d(Old->V)));

        New->f = f_inf(Old->V) + (Old->f - f_inf(Old->V)) * exp(-dt * (alpha_f(Old->V) + beta_f(Old->V)));

        New->x = x_inf(Old->V) + (Old->x - x_inf(Old->V)) * exp(-dt * (alpha_x(Old->V) + beta_x(Old->V)));

        // membrane potential calc                  
        New->V = Old->V + dt / Cm * ( TotalIonCurrent_CPU(0 /* 0 --- index, stands for the single cell model */ , 
                                 Old->V, Old->m, Old->h, Old->J, Old->d, Old->f, Old->x,
                                 (CurrentsOld->INa), (CurrentsOld->IK), (CurrentsOld->IX), 
                                 (CurrentsOld->IS), ConcentrationsOld->Ca)
                                 + IsStim(tCurrent)*I_Stim(0, 0, 2.1) ); // "standart" I_stim = 1e

        // concentrations calc
        ConcentrationsNew->Ca = ConcentrationsOld->Ca + dt * (
                                                            -1e-7*CurrentsOld->IS[0] 
                                                            + 0.07*(1e-7 - ConcentrationsOld->Ca) 
                                                            ); // index "0" --- for array of length 1

        tCurrent += dt;
        //stepNumber += 1;

        // swapping time layers
        struct State* tmp;
        tmp = Old; Old = New; New = tmp;

        struct Concentrations* tmpConc; // for swapping CONCENTRATIONS-type structs
        tmpConc = ConcentrationsOld; ConcentrationsOld = ConcentrationsNew; ConcentrationsNew = tmpConc;


    } // while

    //printf("Second loop is finished; VOfPhase: %.1f mV\n", VOfPhase);
    //std::cin.get();
    // "return" --- is within the loop (look up)
    return stateOfPhase;
}


// we use big "J" for denoting gate var "j": letter "j" is used as a counter in the loops below
void SetInitialConditions_CPU(real* V, real* m, real* h, real* J, real* d, 
real* f, real* x, real* concCa, real value, int numPointsX, int numPointsY) {
    int idx;
    //std::srand(unsigned(1.)); // initial seed for random number generator
    real randomNumber;

    // single initial peak
    //int iCenter = (int)((real)numSegmentsX /2.);
    //int jCenter = (int)((real)numSegmentsY /2.); // tmp values
    //int idxCenter = CalculateLinearCoordinate_CPU(iCenter, jCenter);

    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++) {

            int idxCenter = CalculateLinearCoordinate_CPU(i, j, numPointsX);
            //randomNumber =  ((real)(std::rand() % 20))/20.; // 4phase setting

            // the borders: Dirichlet boundary conditions
            //if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsY - 1)) {
                // TODO: find out about the values
                
                // TODO: wrong formulas below ////////////////////////////////////////////////////
                //int jTilde = numPointsY/2 - j; // j = j^{star}; jTilde = jTilde_{star}
                //int iTilde = i - (numPointsX/2); // i = i^{star}; iTilde = iTilde_{star}

                // for phase calculation: using angle in polar coords
                real LTotal = numPointsX*hx; // should = numPointsY*hy
                real L = (j*hy + hy/2.) - (numPointsY*hy/2.); // LTotal/2. - (j*hx + hx/2.) --- old formula, when j-order was incorrect
                real lsmall = LTotal/2. - ( (numPointsX - 1 - i)*hx + hx/2.);

                real phase = atan2(L, lsmall); // = angle in polar coords; use atan2() func instead of atan() !
                
                // check sign: atan2() returns vals. from [-pi, pi]
                if (phase < 0)
                {
                    phase += 2*M_PI;
                }
                

                //real phaseShifted = phase - M_PI/2.; // phase from R.Syunyaev article
                
                //printf("Phase = %.2f deg.\n", phase*180/M_PI);
                //std::cin.get();
                // TODO //////////////////////////////////////////////////////////////

                // the func returns a std::map of all the vars' values
                ///* std::map<std::string, real> */ real* stateForPhase = VviaPhase(phase);

                //printf("Phase: %.2f deg., VOfPhase = %.2f\n", phase*180./M_PI, stateForPhase["V"]);
                //std::cin.get();

                V[idxCenter] = -60.; //stateForPhase[V_];  //VviaPhase(phase); //M_PI/12. //VRest;
                m[idxCenter] =  0.5; //stateForPhase[m_]; //0.067;//m_inf_CPU(VRest); // 0.5
                h[idxCenter] = 0.5; //stateForPhase[h_]; //0.999; //h_inf_CPU(VRest); // 0.5
                J[idxCenter] = 0.5; //stateForPhase[J_]; //0.;// 0.1; // may be false; TODO perform calcs with higher T 
                d[idxCenter] = 0.5; //stateForPhase[d_]; //0.; //0.;
                f[idxCenter] = 0.5; //stateForPhase[f_]; //1.; //1.;
                x[idxCenter] = 0.5; //stateForPhase[x_];
                concCa[idxCenter] = 3e-7; //stateForPhase[concCa_];

                // for progress checking: in percents
                printf("Set. initial cond: %.2f percent completed\n", 
                        100.*idxCenter / CalculateLinearCoordinate_CPU(numPointsX - 1, numPointsY - 1, numPointsX));
            }

    // after filling the whole area: "fill" borders wiht Neumann boundary cond.
    // the borders: Neumann boundary conditions
    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++)
        {
            int idxCenter = CalculateLinearCoordinate_CPU(i, j, numPointsX);
            
            // borrder cells, including corner cells
            if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsY - 1))
            {
                int idxNear;

                if ((i == 0)) //&& (j >= 1) && (j <= numSegmentsY - 1)) // left border, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i + 1, j, numPointsX);
                if ((j == 0)) //&& (i >= 1) && (i <= numSegmentsX - 1)) // bottom, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i, j + 1, numPointsX);
                if ((j == numPointsY - 1)) // && (i >= 1) && (i <= numSegmentsX - 1)) // top, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i, j - 1, numPointsX);
                if ((i == numPointsX - 1)) // && (j >= 1) && (j <= numSegmentsY - 1)) // right, except for corner cells
                    idxNear = CalculateLinearCoordinate_CPU(i - 1, j, numPointsX);

                // what about corner cells? for now, they are not treated (?)
                V[idxCenter] = V[idxNear];
                m[idxCenter] = m[idxNear];
                h[idxCenter] = h[idxNear];
                J[idxCenter] = J[idxNear];
                d[idxCenter] = d[idxNear];
                f[idxCenter] = f[idxNear];
                x[idxCenter] = x[idxNear];
                concCa[idxCenter] = concCa[idxNear];
            }
        }
}




int main(int argc, char** argv) {

    // setting a GPU for the computations; NOTE: req "openacc.h"!
    //acc_set_device_num(1, acc_device_nvidia);

    // reading the params from the console
    int numSegmentsX = atoi(argv[1]);
    int numSegmentsY = atoi(argv[2]);
    //int serieOfLaunchesNum = atoi(argv[3]);
    
    // storing output file's name in char[]
    /* string */ char tmp_output_file[256];
    sprintf(tmp_output_file, argv[3]); 
    
    int numPointsX = numSegmentsX + 1;
    int numPointsY = numSegmentsY + 1;

    int numPointsTotal = numPointsX * numPointsY;

    // allocating memory
    
    // C++11 style alloc
    //std::vector<real> vecVOld(numPointsTotal);
    //real* VOld = &vecVOld.front();

    // gating vars
    real* VOld = MemAlloc(numPointsTotal); //new real[numPointsTotal];
    real* mOld = MemAlloc(numPointsTotal); 
    real* hOld = MemAlloc(numPointsTotal);
    real* JOld = MemAlloc(numPointsTotal);
    real* dOld = MemAlloc(numPointsTotal);
    real* fOld = MemAlloc(numPointsTotal);
    real* xOld = MemAlloc(numPointsTotal);

    real* VNew = MemAlloc(numPointsTotal);
    real* mNew = MemAlloc(numPointsTotal);
    real* hNew = MemAlloc(numPointsTotal);
    real* JNew = MemAlloc(numPointsTotal);
    real* dNew = MemAlloc(numPointsTotal);
    real* fNew = MemAlloc(numPointsTotal);
    real* xNew = MemAlloc(numPointsTotal);

    // currents
    real* INaOld = MemAlloc(numPointsTotal);
    real* IKOld = MemAlloc(numPointsTotal);
    real* IXOld = MemAlloc(numPointsTotal);
    real* ISOld = MemAlloc(numPointsTotal);

    real* INaNew = MemAlloc(numPointsTotal);
    real* IKNew = MemAlloc(numPointsTotal);
    real* IXNew = MemAlloc(numPointsTotal);
    real* ISNew = MemAlloc(numPointsTotal);

    // concentrations
    real* concCaOld = MemAlloc(numPointsTotal);
    real* concCaNew = MemAlloc(numPointsTotal);

    real* tmp; // a pointer for swapping time-layers 'n' and 'n+1'
    real *tmpConc; // a pointer for swapping time-layers 'n' and 'n+1' for concentrations

    // for output in a loop
    //std::map<std::string, real*> variables;
    real* variables[12]; // number of vars = 12
    variables[V_] = VOld;
    variables[m_] = mOld;
    variables[h_] = hOld;
    variables[J_] = JOld;
    variables[d_] = dOld;
    variables[f_] = fOld;
    variables[x_] = xOld;
    
    variables[INa_] = INaOld;
    variables[IK_] = IKOld;
    variables[IX_] = IXOld;
    variables[IS_] = ISOld;

    variables[concCa_] = concCaOld;

    // = {"V", "m", "h", "J", "d", "f", "x"};


    // initializing before timesteppin'
    SetInitialConditions_CPU(VOld, mOld, hOld, JOld, dOld, fOld, xOld, concCaOld, 0., numPointsX, numPointsY);
    //SetInitialConditions_CPU(VNew, mNew, hNew, JNew, dNew, fNew, xNew, 0.); // for avoiding "junk" values in all '...New' arrays

    real tCurrent = 0.;
    int stepNumber = 0;
    int counterOutput = 1;

    printf("Timesteppin' begins...\n");
    clock_t start = clock();

// pragmas without "-acc" flag --- are ignored?
#pragma acc data copy(VOld[0:numPointsTotal], mOld[0:numPointsTotal], hOld[0:numPointsTotal], \
JOld[0:numPointsTotal], dOld[0:numPointsTotal], fOld[0:numPointsTotal], xOld[0:numPointsTotal], \
INaOld[0:numPointsTotal], IKOld[0:numPointsTotal], IXOld[0:numPointsTotal], ISOld[0:numPointsTotal], \
concCaOld[0:numPointsTotal], \
VNew[0:numPointsTotal], mNew[0:numPointsTotal], hNew[0:numPointsTotal], \
JNew[0:numPointsTotal], dNew[0:numPointsTotal], fNew[0:numPointsTotal], xNew[0:numPointsTotal], \
concCaNew[0:numPointsTotal]) \
deviceptr(tmp, tmpConc)
{
    // main loop: timestepping
    while (tCurrent < T) 
    {

        // TODO: change order of indexing (i, j)
    
    // dont forget to try "kernels" keyword
	#pragma acc parallel \
	present(VOld[0:numPointsTotal], mOld[0:numPointsTotal], hOld[0:numPointsTotal], \
    JOld[0:numPointsTotal], dOld[0:numPointsTotal], fOld[0:numPointsTotal], xOld[0:numPointsTotal], \
    concCaOld[0:numPointsTotal], \
    VNew[0:numPointsTotal], mNew[0:numPointsTotal], hNew[0:numPointsTotal], \
    JNew[0:numPointsTotal], dNew[0:numPointsTotal], fNew[0:numPointsTotal], xNew[0:numPointsTotal], \
    concCaNew[0:numPointsTotal])
	{
	
	#pragma acc loop collapse(2) independent
	for (int j = 0; j < numPointsY; j++)
            for (int i = 0; i < numPointsX; i++) 
            {

                int idxCenter = CalculateLinearCoordinate(i, j, numPointsX);
                
                // inner cells
                if (i >= 1 && j >= 1 && i <= (numSegmentsX - 1) && j <= (numSegmentsY - 1))
                {
                    // for short names
                    int idxUp = CalculateLinearCoordinate(i, j + 1, numPointsX);
                    int idxDown = CalculateLinearCoordinate(i, j - 1, numPointsX);
                    int idxLeft = CalculateLinearCoordinate(i - 1, j, numPointsX);
                    int idxRight = CalculateLinearCoordinate(i + 1, j, numPointsX);

                    
                    ///////////////// gating variables: ode ("reaction") step

                    // TODO: make only ONE read of VOld[idxCenter], etc from memory; to more speedup, esp. for GPU
                    mNew[idxCenter] = m_inf(VOld[idxCenter]) + (mOld[idxCenter] - m_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_m(VOld[idxCenter]) + beta_m(VOld[idxCenter])));

                    
                    hNew[idxCenter] = h_inf(VOld[idxCenter]) + (hOld[idxCenter] - h_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_h(VOld[idxCenter]) + beta_h(VOld[idxCenter])));
                    
                    JNew[idxCenter] = j_inf(VOld[idxCenter]) + (JOld[idxCenter] - j_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_j(VOld[idxCenter]) + beta_j(VOld[idxCenter])));
                   

                    dNew[idxCenter] = d_inf(VOld[idxCenter]) + (dOld[idxCenter] - d_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_d(VOld[idxCenter]) + beta_d(VOld[idxCenter])));
                    
                    fNew[idxCenter] = f_inf(VOld[idxCenter]) + (fOld[idxCenter] - f_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_f(VOld[idxCenter]) + beta_f(VOld[idxCenter])));
                   
                    xNew[idxCenter] = x_inf(VOld[idxCenter]) + (xOld[idxCenter] - x_inf(VOld[idxCenter]))
                                                                * exp(-dt * (alpha_x(VOld[idxCenter]) + beta_x(VOld[idxCenter])));
                    
                    
                    /*
                    // for steady state's calculation
                    VNew[idxCenter] = VRest;
                    */

                    
                    //////////////////
                    // "discrete diffusion" step
                    VNew[idxCenter] = VOld[idxCenter]
                    // uncomment if cells are connected; otherwize --- Nans
                     + dt / Cm * (
                            Dx  * (VOld[idxRight] - 2 * VOld[idxCenter] + VOld[idxLeft])
                            + Dy  * (VOld[idxUp] - 2 * VOld[idxCenter] + VOld[idxDown])
                                )
                    + dt / Cm * (
                                                        TotalIonCurrent(idxCenter, VOld[idxCenter], mOld[idxCenter],
                                                             hOld[idxCenter], JOld[idxCenter], 
                                                             dOld[idxCenter], fOld[idxCenter], xOld[idxCenter],
                                                            INaOld, IKOld, IXOld, ISOld, concCaOld[idxCenter])
                                                                        /*+ 0.*IsStim(tCurrent)*I_Stim(i, j, 2.1)*/
                                                ); // "standart" I_stim = 1e0;
                    
                    
                    
                    // reaction step
                    /* VNew[idxCenter] +=  dt / Cm * (
                                                        TotalIonCurrent(idxCenter, VOld[idxCenter], mOld[idxCenter],
                                                             hOld[idxCenter], JOld[idxCenter], 
                                                             dOld[idxCenter], fOld[idxCenter], xOld[idxCenter],
                                                            INaOld, IKOld, IXOld, ISOld, concCaOld[idxCenter])
                                                                        + IsStim(tCurrent)*I_Stim(i, j, 5.)
                                                ); // "standart" I_stim = 1e0;
                    */
                    
                    // concentrations calc
                    concCaNew[idxCenter] = concCaOld[idxCenter] + dt * (
                                                            -1.*1e-7*ISOld[idxCenter] 
                                                            + 0.07*(1e-7 - concCaOld[idxCenter]) 
                                                            ); // index "0" --- for array of length 1
                    
               } // if
               
               // the borders: Neumann boundary conditions
               else
               {
                    int idxNear;
                    
                    if ((i == 0) && (j >= 1) && (j <= numSegmentsY - 1)) // left border, except for corner cells
                        idxNear = CalculateLinearCoordinate(i + 1, j, numPointsX);
                    else if ((j == 0) && (i >= 1) && (i <= numSegmentsX - 1)) // bottom, except for corner cells
                        idxNear = CalculateLinearCoordinate(i, j + 1, numPointsX);
                    else if ((j == numSegmentsY) && (i >= 1) && (i <= numSegmentsX - 1)) // top, except for corner cells
                        idxNear = CalculateLinearCoordinate(i, j - 1, numPointsX);
                    else if ((i == numSegmentsX) && (j >= 1) && (j <= numSegmentsY - 1)) // right, except for corner cells
                        idxNear = CalculateLinearCoordinate(i - 1, j, numPointsX);
                    else { // if corner cell
                        continue; // do nothing, continue the "i,j" loop
                    }

                    // what about corner cells? for now, they are not treated (?)
                    // Neumann boundary cond setting
                    VNew[idxCenter] = VNew[idxNear];
                    mNew[idxCenter] = mNew[idxNear];
                    hNew[idxCenter] = hNew[idxNear];
                    JNew[idxCenter] = JNew[idxNear];
                    dNew[idxCenter] = dNew[idxNear];
                    fNew[idxCenter] = fNew[idxNear];
                    xNew[idxCenter] = xNew[idxNear];
                    concCaNew[idxCenter] = concCaNew[idxNear];
               }

            } // for
	
	} // acc kernels
    
    if ((stepNumber % 4000) == 0) // output each 10 msec: 10/dt(=0.005 ms) = 2000 (old val.)
    { // output each 10 msec: 10/dt = 2000 (old val.)
        //if ( (stepNumber) % (int)(T/dt/500)  == 0 ) {
        #pragma acc update host(VOld[0:numPointsTotal])
            variables[V_] = VOld;
            //variables[m_] = mOld;
            //variables[h_] = hOld;
            //variables[J_] = JOld;
            //variables[d_] = dOld;
            //variables[f_] = fOld;
            //variables[x_] = xOld;
            //variables[INa_] = INaOld;
            //variables[IK_] = IKOld;
            //variables[IX_] = IXOld;
            //variables[IS_] = ISOld;
            //variables[concCa_] = concCaOld;

            int outNumber = stepNumber;
                
            //if (variable.first.compare("V") == 0 ) // output only "V"
            //{
            Write2VTK(numPointsX, variables[V_], hx, outNumber); // for now: numPointsX == numPointsY
            //Write2VTK("V", numPointsX, variables["V"], hx, counterOutput); // for now: numPointsX == numPointsY
            //}
            //}
            printf("Time: %.2f percent completed\n", 100.*stepNumber*dt/T);
	        counterOutput++;
        }
        
        tCurrent += dt;
        stepNumber += 1;

        // swapping time-layers
        ////// swap V
        tmp = VOld; VOld = VNew; VNew = tmp;
        ///// swap m
        tmp = mOld; mOld = mNew; mNew = tmp;
        ///// swap h
        tmp = hOld; hOld = hNew; hNew = tmp;
        ///// swap d
        tmp = JOld; JOld = JNew; JNew = tmp;
        ///// swap d
        tmp = dOld; dOld = dNew; dNew = tmp;
        ///// swap f
        tmp = fOld; fOld = fNew; fNew = tmp;
        ///// swap x
        tmp = xOld; xOld = xNew; xNew = tmp;

        ////// swap concCa
        tmpConc = concCaOld; concCaOld = concCaNew; concCaNew = tmpConc;

    } // while (tCurrent < T)


} // acc data

    real elapsedTime = (real)( ((real)(clock() - start))/CLOCKS_PER_SEC );
    printf("\nCalculations finished. Elapsed time = %.2e sec\n", elapsedTime);

    // printing elapsed time into a file
    FILE* ff = fopen(tmp_output_file, "w");
    fprintf(ff, "%.2f", elapsedTime);
    fclose(ff);
    
    // cleaning up
    free(VOld);
    free(VNew);
    free(mOld);
    free(mNew);
    free(hOld);
    free(hNew);
    
    free(JOld);
    free(JNew);
    
    free(dOld);
    free(dNew);
    free(fOld);
    free(fNew);

    free(xOld);
    free(xNew);

    free(concCaOld);
    free(concCaNew);


    return 0;
}
