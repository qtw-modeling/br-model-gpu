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
#include "extra.h"
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
#define T (25.) //(1000.) // old val: 500 // endtime (in ms)
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
real* CalculateStateFromFile(real phase) 
{
    
    FILE* fp;
    /* открытие для чтения */
    if( (fp = fopen("phase_data_br_model.bin", "rb")) == NULL) 
    {
        printf("Cannot open file");
        exit(1);
    }

    real phaseLocalOld, phaseLocalNew;
    real* stateOfPhaseOld = (real*)malloc(8*sizeof(real));
    real* stateOfPhaseNew = (real*)malloc(8*sizeof(real));
    real* stateOfPhaseFinal = (real*)malloc(8*sizeof(real));
    real* tmp; // 4switching "timelayers"

    // setting "initial conditions"
    fread(&phaseLocalOld, sizeof(real), 1, fp);
    /* reading the whole array in one step */
    fread(stateOfPhaseOld, 8*sizeof(real), 1, fp);
    
    // going over all data in the bin file
    while (1)
    {
        
        fread(&phaseLocalNew, sizeof(real), 1, fp);
        /* reading the whole array in one step */
        fread(stateOfPhaseNew, 8*sizeof(real), 1, fp);
        
        if (phaseLocalNew >= phase)
        {
            //break; // stands for right-neighbour interpolation, since we have "passed" over "phase"
            
            // using linear interpolation
            for (int ii = 0; ii <= 7; ii++)
                stateOfPhaseFinal[ii] = CalculateLinearInterpolate(phase, phaseLocalOld, phaseLocalNew, stateOfPhaseOld[ii], stateOfPhaseNew[ii]);
        
            break;
        }

        tmp = stateOfPhaseOld; stateOfPhaseOld = stateOfPhaseNew; stateOfPhaseNew = tmp;
        

    }
    
    fclose(fp);
    
    return stateOfPhaseFinal;
}


void SetInitialConditions_CPU(real* V, real* m, real* h, real* J, real* d, real* f, 
real* x, real* concCa, real value, int numPointsX, int numPointsY) 
{
    int idx;
    //std::srand(unsigned(1.)); // initial seed for random number generator
    //real randomNumber;

    // single initial peak
    //int iCenter = (int)((real)numSegmentsX /2.);
    //int jCenter = (int)((real)numSegmentsY /2.); // tmp values
    //int idxCenter = CalculateLinearCoordinate_CPU(iCenter, jCenter);

    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++) 
        {

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

                real* stateForPhase = CalculateStateFromFile(phase);

                //printf("Phase: %.2f deg., VOfPhase = %.2f\n", phase*180./M_PI, stateForPhase["V"]);
                //std::cin.get();

                V[idxCenter] = stateForPhase[V_];  //VviaPhase(phase); //M_PI/12. //VRest;
                m[idxCenter] = stateForPhase[m_]; //0.067;//m_inf_CPU(VRest); // 0.5
                h[idxCenter] = stateForPhase[h_]; //0.999; //h_inf_CPU(VRest); // 0.5
                J[idxCenter] = stateForPhase[J_]; //0.;// 0.1; // may be false; TODO perform calcs with higher T 
                d[idxCenter] = stateForPhase[d_]; //0.; // may be false; TODO perform calcs with higher T  
                f[idxCenter] = stateForPhase[f_]; //0.; //0.;
                x[idxCenter] = stateForPhase[x_]; //1.; //1.;
                concCa[idxCenter] = stateForPhase[concCa_];

                // for progress checking: in percents
                //printf("Set. initial cond: %.2f percent completed\n", 
                //        100.*idxCenter / CalculateLinearCoordinate_CPU(numPointsX - 1, numPointsY - 1, numPointsX));
            }


    // after filling the whole area: "fill" borders wiht Neumann boundary cond.
    // the borders: Neumann boundary conditions
    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++)
        {
            int idxCenter = CalculateLinearCoordinate_CPU(i, j, numPointsX);
            
            // borrder cells, including corner cells
            if (i == 0 || j == 0 || i == (numPointsX - 1) || j == (numPointsX - 1))
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

    // we pass number of cells as command line args; 
    // reading the params from the console:
    int numSegmentsX = atoi(argv[1]) + 1;
    int numSegmentsY = atoi(argv[2]) + 1;
    //int serieOfLaunchesNum = atoi(argv[3]);
    
    // storing output file's name in char[]
    /* string */ char tmp_output_file[256];
    sprintf(tmp_output_file, argv[3]); 
    
    int numPointsX = numSegmentsX + 1; // includes 2 ghost cells; numCells + 2 = numPointsX
    int numPointsY = numSegmentsY + 1; // same

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
    //SetInitialConditions_CPU(VNew, mNew, hNew, JNew, dNew, fNew, xNew, concCaOld, 0., numPointsX, numPointsY); // for avoiding "junk" values in all '...New' arrays

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
    concCaNew[0:numPointsTotal]) \
    vector_length(32) num_workers(1) //num_gangs(32)
	{
	
    // loop over ALL cells
	//#pragma acc loop independent vector(32) worker(2) gang(256)
    #pragma acc loop gang // as many gangs (= blocks) as needed
    for (int j = 0; j < numPointsY; j++)
    {       
        //#pragma acc loop independent vector(32) worker(2) gang(256)
        #pragma acc loop vector
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
        } // for i
    } // for j

                
               // the borders: Neumann boundary conditions
                #pragma acc loop gang
	            for (int j = 0; j < numPointsY; j++)
                {
                    #pragma acc loop vector
                    for (int i = 0; i < numPointsX; i++) 
                    {
                        int idxCenter = CalculateLinearCoordinate(i, j, numPointsX);
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
                    } // for i
                } // for j
	
	} // acc kernels
    
    if ((stepNumber % 2000) == 0) // output each 10 msec: 10/dt(=0.005 ms) = 2000 (old val.)
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
