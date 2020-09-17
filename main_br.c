//
// Created by QuoTheWhite on 27/03/2019.
//#include "openacc.h"
// 4 C++11 includes
//#include <vector>

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
//#define T (25.) //(1000.) // old val: 500 // endtime (in ms)
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



int CalculateLinearCoordinate_CPU(int i, int j, int numPointsX) {
    return i + j*numPointsX;
}

#pragma acc routine
int CalculateLinearCoordinate(int i, int j, int numPointsX) {
    return i + j*numPointsX;
}

#pragma acc routine seq
void /* real* */ SliceGhosts2D(int xDimWithGhosts, int yDimWithGhosts, real* arrWithGhosts, real* arrNoGhosts)
{
    //real* arrNoGhosts = MemAlloc( (xDimWithGhosts - 2)*(yDimWithGhosts - 2) ); // exclude 2 ghost cells in each DIM
    
    // aux. vars
    int idx, counter, J, offset;
    counter = 0;
    offset = xDimWithGhosts;
    
    // filling arrNoGhosts consequently 
    //#pragma acc loop vector
    /*
    for (int k = 0; k < (xDimWithGhosts - 2)*(yDimWithGhosts - 2); k++)
    {
        J = 2*counter + 1; // just a formula for an odd num
        idx = offset + k + J; 
        arrNoGhosts[k] = arrWithGhosts[idx];
        
        if ( (k + 1) % (xDimWithGhosts - 2) == 0) // TODO: check correctness of indexing of arraWithGhosts[idx] by hand!
            counter += 1;
    }
    */

    int k = 0; // for indexing arrNoGhosts (its 1D-arrray)
    while (k < xDimWithGhosts - 2) // without 2 ghost cells
    {
        // output without ghost cells
        for (int j = 1; j < xDimWithGhosts - 1; j++) 
        {
            for (int i = 1; i < xDimWithGhosts - 1; i++)
            {
                arrNoGhosts[k] = arrWithGhosts[j * xDimWithGhosts + i]; // [linear index; should be same as in CalcLinearCoord(...)]
                k += 1; // moving to next element of arrNoGhosts
            }
        }
    }

    //return arrNoGhosts;
}


real FuncSinglePeriod(real t)
{
    // We consider onlt t > 0! //

    if (t < PeriodStim / 2.)
        return 0.; // first: no stimulation
    else
        return 1.; // second: stimulation


/*
 // "__П__ " --- func looks like this
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

    for (int j = 0; j < numPointsY; j++)
        for (int i = 0; i < numPointsX; i++) 
        {

            int idxCenter = CalculateLinearCoordinate_CPU(i, j, numPointsX);

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
                // TODO //////////////////////////////////////////////////////////////

                real* stateForPhase = CalculateStateFromFile(phase);

                //printf("Phase: %.2f deg., VOfPhase = %.2f\n", phase*180./M_PI, stateForPhase["V"]);
                

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
    int numCellsX = atoi(argv[1]);
    int numCellsY = atoi(argv[2]);
    int numCellsTotal = numCellsX*numCellsY;
    //int serieOfLaunchesNum = atoi(argv[3]);
    
    const int T = atoi(argv[3]);

    // storing output file's name in char[]
    /* string */ char tmp_output_file[256];
    sprintf(tmp_output_file, argv[4]);
    
    int numPointsX = numCellsX + 2; // includes 2 ghost cells; numCells + 2 = numPointsX
    int numPointsY = numCellsY + 2; // same
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

    real* printV = MemAlloc(numCellsTotal);

    // for output in a loop
    //std::map<std::string, real*> variables;
    const int NUM_VARS = 12;
    
    real* variablesOld[NUM_VARS]; // number of vars = 12
    variablesOld[V_] = VOld;
    variablesOld[m_] = mOld;
    variablesOld[h_] = hOld;
    variablesOld[J_] = JOld;
    variablesOld[d_] = dOld;
    variablesOld[f_] = fOld;
    variablesOld[x_] = xOld;
    
    variablesOld[INa_] = INaOld;
    variablesOld[IK_] = IKOld;
    variablesOld[IX_] = IXOld;
    variablesOld[IS_] = ISOld;

    variablesOld[concCa_] = concCaOld;

    // array of vars (4 deletion within a loop in the end)
    real* variablesNew[NUM_VARS]; // number of vars = 12
    variablesNew[V_] = VNew;
    variablesNew[m_] = mNew;
    variablesNew[h_] = hNew;
    variablesNew[J_] = JNew;
    variablesNew[d_] = dNew;
    variablesNew[f_] = fNew;
    variablesNew[x_] = xNew;
    
    variablesNew[INa_] = INaNew;
    variablesNew[IK_] = IKNew;
    variablesNew[IX_] = IXNew;
    variablesNew[IS_] = ISNew;

    variablesNew[concCa_] = concCaNew;

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
create(printV[0:numCellsTotal]) \
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
	
    // loop over inner cells
	//#pragma acc loop independent vector(32) worker(2) gang(256)
    #pragma acc loop gang // as many gangs (= blocks) as needed
    for (int j = 1; j <= numPointsY - 1; j++)
    {       
        //#pragma acc loop independent vector(32) worker(2) gang(256)
        #pragma acc loop vector
        for (int i = 1; i <= numPointsX - 1; i++)
        {

            int idxCenter = CalculateLinearCoordinate(i, j, numPointsX);
                
            // inner cells
            //if (i >= 1 && j >= 1 && i <= (numSegmentsX - 1) && j <= (numSegmentsY - 1))
            //{
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
                    
            //   } // if
        } // for i
    } // for j

    } // acc parallel 4 inner cells

            
            // setting Neumann BCs, NEW VERS;
            // corner cells are not treated. 
                //#pragma acc parallel async(0) \
                //present(VNew[0:numPointsTotal])
                //{
                    #pragma acc parallel async(0) \
                    present(VNew[0:numPointsTotal])
                    {
                    #pragma acc loop vector // seq
                    for (int j = 1; j < numPointsY - 1; j++)
                    {
                        int idxCenterLeftBord = CalculateLinearCoordinate(0, j, numPointsX);
                        int idxCenterRightBord = CalculateLinearCoordinate(numPointsX - 1, j, numPointsX);
                        
                        int idxNearLeft = CalculateLinearCoordinate(1, j, numPointsX);
                        int idxNearRight = CalculateLinearCoordinate(numPointsX - 1 - 1, j, numPointsX);

                        VNew[idxCenterLeftBord] = VNew[idxNearLeft];
                        VNew[idxCenterRightBord] = VNew[idxNearRight];
                    }
                    }

                    #pragma acc parallel async(1) \
                    present(VNew[0:numPointsTotal])
                    {
                    #pragma acc loop vector // seq
                    for (int i = 1; i < numPointsX - 1; i++)
                    {
                        int idxCenterBottomBord = CalculateLinearCoordinate(i, 0, numPointsX);
                        int idxCenterTopBord = CalculateLinearCoordinate(i, numPointsY - 1, numPointsX);
                        
                        int idxNearBottom = CalculateLinearCoordinate(i, 1, numPointsX);
                        int idxNearTop = CalculateLinearCoordinate(i, numPointsY - 1 - 1, numPointsX);

                        VNew[idxCenterBottomBord] = VNew[idxNearBottom];
                        VNew[idxCenterTopBord] = VNew[idxNearTop];
                    }
                    }
                
                #pragma acc wait(0, 1)



                /*
               // the borders: Neumann boundary conditions, OLD VERS
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
                        // Neumann boundary cond setting;
                        // we need only set BC for V, not other vars?
                        VNew[idxCenter] = VNew[idxNear];
                        //mNew[idxCenter] = mNew[idxNear];
                        //hNew[idxCenter] = hNew[idxNear];
                        //JNew[idxCenter] = JNew[idxNear];
                        //dNew[idxCenter] = dNew[idxNear];
                        //fNew[idxCenter] = fNew[idxNear];
                        //xNew[idxCenter] = xNew[idxNear];
                        //concCaNew[idxCenter] = concCaNew[idxNear];
                    } // for i
                } // for j
                */

	//} // acc parallel 4 inner cells
    
    if ((stepNumber % 2000) == 0) // output each 10 msec: 10/dt(=0.005 ms) = 2000 (old val.)
    { 
        #pragma acc parallel \
        present(VOld[0:numPointsTotal], printV[0:numCellsTotal])
        {
            SliceGhosts2D(numPointsX, numPointsY, VOld, printV);
        }
        
        // output each 10 msec: 10/dt = 2000 (old val.)
        //if ( (stepNumber) % (int)(T/dt/500)  == 0 ) {
        #pragma acc update host(printV[0:numCellsTotal]) // (VOld[0:numPointsTotal])
            
            //variablesOld[V_] = printV;
            
            //variablesOld[V_] = VOld;
            //variablesOld[m_] = mOld;
            //variablesOld[h_] = hOld;
            //variablesOld[J_] = JOld;
            //variablesOld[d_] = dOld;
            //variablesOld[f_] = fOld;
            //variablesOld[x_] = xOld;
            //variablesOld[INa_] = INaOld;
            //variablesOld[IK_] = IKOld;
            //variablesOld[IX_] = IXOld;
            //variablesOld[IS_] = ISOld;
            //variablesOld[concCa_] = concCaOld;

            //real* printV = SliceGhosts2D(numPointsX, numPointsY, variablesOld[V_]);

            int outNumber = stepNumber;
                
            //if (variable.first.compare("V") == 0 ) // output only "V"
            //{
            
            Write2VTK(numCellsX, printV, hx, outNumber);
            //Write2VTKWithGhosts(numPointsX, variablesOld[V_], hx, outNumber); // for now: numPointsX == numPointsY
            //Write2VTK_noGhosts(numPointsX, variablesOld[V_], hx, outNumber); // for now: numPointsX == numPointsY
            //Write2VTK("V", numPointsX, variablesOld["V"], hx, counterOutput); // for now: numPointsX == numPointsY
            //}
            //}
            printf("Progress: %.2f %% completed;\n", 100.*stepNumber*dt/T);
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
    printf("\nCalculations finished. Elapsed time = %.2e sec.\n", elapsedTime);

    // printing elapsed time into a file
    FILE* ff = fopen(tmp_output_file, "w");
    fprintf(ff, "%.2f", elapsedTime);
    fclose(ff);
    
    // cleaning up
    for (int i = 0; i < NUM_VARS; i++)
    {
        free(variablesOld[i]);
        free(variablesNew[i]);
    }


    return 0;
}
