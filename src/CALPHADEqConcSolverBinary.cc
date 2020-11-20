#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADFunctions.h"

#include <cassert>
#include <cmath>
#include <iostream>

//=======================================================================

void CALPHADEqConcentrationSolverBinary::RHS(
    const double* const c, // composition of species A in various phases
    double* const fvec)
{
    // tbox::pout<<"Compute RHS for CALPHAD..."<<endl;
    double dfdci[2];

    // loop over phases
    for (int ii = 0; ii < 2; ii++)
    {
        dfdci[ii] = fA_[ii] - fB_[ii]
                    + CALPHADcomputeFMix_derivBinary(
                          L0_[ii], L1_[ii], L2_[ii], L3_[ii], c[ii])
                    + CALPHADcomputeFIdealMix_derivBinary(RT_, c[ii]);
    }

    fvec[0] = c[0] * fA_[0] + (1.0 - c[0]) * fB_[0]
              + CALPHADcomputeFIdealMixBinary(RT_, c[0])
              + CALPHADcomputeFMixBinary(L0_[0], L1_[0], L2_[0], L3_[0], c[0])
              - c[1] * fA_[1] - (1.0 - c[1]) * fB_[1]
              - CALPHADcomputeFIdealMixBinary(RT_, c[1])
              - CALPHADcomputeFMixBinary(L0_[1], L1_[1], L2_[1], L3_[1], c[1])
              - (c[0] - c[1]) * (dfdci[1]);

    // dfL/dcL - dfS/dcS
    fvec[1] = dfdci[0] - dfdci[1];
}

//=======================================================================

void CALPHADEqConcentrationSolverBinary::Jacobian(
    const double* const c, double** const fjac)
{
    // tbox::pout<<"Compute Jacobian for CALPHAD..."<<endl;
    double dfdci[2];

    // loop over phases
    for (int ii = 0; ii < 2; ii++)
    {

        dfdci[ii] = fA_[ii] - fB_[ii]
                    + CALPHADcomputeFMix_derivBinary(
                          L0_[ii], L1_[ii], L2_[ii], L3_[ii], c[ii])
                    + CALPHADcomputeFIdealMix_derivBinary(RT_, c[ii]);
    }

    // f[i][j]=df[i]/dc[j]
    fjac[0][0] = dfdci[0] - dfdci[1];
    fjac[0][1] = -(c[0] - c[1])
                 * (CALPHADcomputeFIdealMix_deriv2Binary(RT_, c[1])
                       + CALPHADcomputeFMix_deriv2Binary(
                             L0_[1], L1_[1], L2_[1], L3_[1], c[1]));

    fjac[1][0]
        = CALPHADcomputeFMix_deriv2Binary(L0_[0], L1_[0], L2_[0], L3_[0], c[0])
          + CALPHADcomputeFIdealMix_deriv2Binary(RT_, c[0]);

    fjac[1][1]
        = -CALPHADcomputeFMix_deriv2Binary(L0_[1], L1_[1], L2_[1], L3_[1], c[1])
          - CALPHADcomputeFIdealMix_deriv2Binary(RT_, c[1]);
}

//=======================================================================

int CALPHADEqConcentrationSolverBinary::ComputeConcentration(double* const conc,
    const double RTinv, const double* const L0, const double* const L1,
    const double* const L2, const double* const L3, const double* const fA,
    const double* const fB)
{
    RTinv_ = RTinv;
    RT_    = 1. / RTinv;

    for (int ii = 0; ii < 2; ii++)
    {
        L0_[ii] = L0[ii];
        L1_[ii] = L1[ii];
        L2_[ii] = L2[ii];
        L3_[ii] = L3[ii];
        fA_[ii] = fA[ii];
        fB_[ii] = fB[ii];
    }

    return DampedNewtonSolver::ComputeSolution(conc, 2);
}
