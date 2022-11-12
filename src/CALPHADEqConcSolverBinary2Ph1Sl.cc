#include "CALPHADEqConcSolverBinary2Ph1Sl.h"
#include "CALPHADFunctions.h"

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
void CALPHADEqConcSolverBinary2Ph1Sl::RHS(
    const double* const c, // composition of species A in two phases
    double* const fvec)
{
    // Get the sublattice site fractions, assuming (A)_p (A,B)_q
    double ypp_A[2];
    for (int i = 0; i < 2; ++i)
    {
        ypp_A[i] = (p_[i] + q_[i]) * c[i] - p_[i];
    }

    double dfdci[2];

    dfdci[0] = fA_[0] - fB_[0]
               + CALPHADcomputeFMix_derivBinary(
                     Lmix_L_[0], Lmix_L_[1], Lmix_L_[2], Lmix_L_[3], ypp_A[0])
               + CALPHADcomputeFIdealMix_derivBinary(q_[0] * RT_, ypp_A[0]);
    dfdci[1] = fA_[1] - fB_[1]
               + CALPHADcomputeFMix_derivBinary(
                     Lmix_A_[0], Lmix_A_[1], Lmix_A_[2], Lmix_A_[3], ypp_A[1])
               + CALPHADcomputeFIdealMix_derivBinary(q_[1] * RT_, ypp_A[1]);

    fvec[0] = (ypp_A[0] * fA_[0] + (1.0 - ypp_A[0]) * fB_[0]
                  + CALPHADcomputeFIdealMixBinary(q_[0] * RT_, ypp_A[0])
                  + CALPHADcomputeFMixBinary(Lmix_L_[0], Lmix_L_[1], Lmix_L_[2],
                        Lmix_L_[3], ypp_A[0]))
                  / (p_[0] + q_[0])
              - (ypp_A[1] * fA_[1] + (1.0 - ypp_A[1]) * fB_[1]
                    + CALPHADcomputeFIdealMixBinary(q_[1] * RT_, ypp_A[1])
                    + CALPHADcomputeFMixBinary(Lmix_A_[0], Lmix_A_[1],
                          Lmix_A_[2], Lmix_A_[3], ypp_A[1]))
                    / (p_[1] + q_[1])
              - (c[0] - c[1]) * dfdci[1];

    // dfL/dcL - dfS/dcS
    fvec[1] = dfdci[0] - dfdci[1];
}

//=======================================================================

void CALPHADEqConcSolverBinary2Ph1Sl::Jacobian(
    const double* const c, JacobianDataType** const fjac)
{
    // Get the sublattice site fractions, assuming (A)_p (A,B)_q
    double ypp_A[2];
    for (int i = 0; i < 2; ++i)
    {
        ypp_A[i] = (p_[i] + q_[i]) * c[i] - p_[i];
    }

    double dfdci[2];

    // loop over phases
    dfdci[0] = fA_[0] - fB_[0]
               + CALPHADcomputeFMix_derivBinary(
                     Lmix_L_[0], Lmix_L_[1], Lmix_L_[2], Lmix_L_[3], ypp_A[0])
               + CALPHADcomputeFIdealMix_derivBinary(q_[0] * RT_, ypp_A[0]);
    dfdci[1] = fA_[1] - fB_[1]
               + CALPHADcomputeFMix_derivBinary(
                     Lmix_A_[0], Lmix_A_[1], Lmix_A_[2], Lmix_A_[3], ypp_A[1])
               + CALPHADcomputeFIdealMix_derivBinary(q_[1] * RT_, ypp_A[1]);

    // f[i][j]=drhS[i]/dc[j]
    fjac[0][0] = dfdci[0] - dfdci[1];
    fjac[0][1] = -(c[0] - c[1]) * (p_[1] + q_[1])
                 * (CALPHADcomputeFIdealMix_deriv2Binary(q_[1] * RT_, ypp_A[1])
                       + CALPHADcomputeFMix_deriv2Binary(Lmix_A_[0], Lmix_A_[1],
                             Lmix_A_[2], Lmix_A_[3], ypp_A[1]));

    fjac[1][0]
        = (p_[0] + q_[0])
          * (CALPHADcomputeFMix_deriv2Binary(
                 Lmix_L_[0], Lmix_L_[1], Lmix_L_[2], Lmix_L_[3], ypp_A[0])
                + CALPHADcomputeFIdealMix_deriv2Binary(q_[0] * RT_, ypp_A[0]));

    fjac[1][1]
        = -(p_[1] + q_[1])
          * (CALPHADcomputeFMix_deriv2Binary(
                 Lmix_A_[0], Lmix_A_[1], Lmix_A_[2], Lmix_A_[3], ypp_A[1])
                + CALPHADcomputeFIdealMix_deriv2Binary(q_[1] * RT_, ypp_A[1]));
}

//=======================================================================

void CALPHADEqConcSolverBinary2Ph1Sl::setup(const double RTinv,
    const CalphadDataType* const Lmix_L, const CalphadDataType* const Lmix_A,
    const CalphadDataType* const fA, const CalphadDataType* const fB,
    const int* const p, const int* const q)
{
    RTinv_ = RTinv;
    RT_    = 1. / RTinv;

    for (int ii = 0; ii < 4; ii++)
    {
        Lmix_L_[ii] = Lmix_L[ii];
        Lmix_A_[ii] = Lmix_A[ii];
    }
    for (int ii = 0; ii < 2; ii++)
    {
        fA_[ii] = fA[ii];
        fB_[ii] = fB[ii];
    }
    for (int ii = 0; ii < 2; ii++)
    {
        p_[ii] = static_cast<double>(p[ii]);
        q_[ii] = static_cast<double>(q[ii]);
    }
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
