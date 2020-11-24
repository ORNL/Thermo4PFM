#include "CALPHADEqConcSolverTernary.h"
#include "CALPHADFunctions.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace Thermo4PFM
{

CALPHADEqConcentrationSolverTernary::CALPHADEqConcentrationSolverTernary()
    : DampedNewtonSolver(4)
{
    fA_[0] = std::nan("");
    fA_[1] = std::nan("");
    fB_[0] = std::nan("");
    fB_[1] = std::nan("");
    fC_[0] = std::nan("");
    fC_[1] = std::nan("");
}

void CALPHADEqConcentrationSolverTernary::RHS(
    const double* const c, double* const fvec)
{
    assert(fA_[0] == fA_[0]);
    assert(fC_[1] == fC_[1]);

    const double* const cL = &c[0]; // composition of Species A and B in phase L
    const double* const cS = &c[2]; // composition of Species A and B in phase S
    // tbox::pout<<"Compute RHS for CALPHAD..."<<endl;

    double derivIdealMixL[2];
    CALPHADcomputeFIdealMix_derivTernary(RT_, cL[0], cL[1], derivIdealMixL);

    double derivFMixL[2];
    CALPHADcomputeFMix_derivTernary(
        L_AB_L_, L_AC_L_, L_BC_L_, L_ABC_L_, cL[0], cL[1], derivFMixL);

    double dfLdciL[2];
    // 1st species
    dfLdciL[0] = fA_[0] - fC_[0] + derivFMixL[0] + derivIdealMixL[0];

    // 2nd species
    dfLdciL[1] = fB_[0] - fC_[0] + derivFMixL[1] + derivIdealMixL[1];

    double derivIdealMixS[2];
    CALPHADcomputeFIdealMix_derivTernary(RT_, cS[0], cS[1], derivIdealMixS);

    double derivFMixS[2];
    CALPHADcomputeFMix_derivTernary(
        L_AB_S_, L_AC_S_, L_BC_S_, L_ABC_S_, cS[0], cS[1], derivFMixS);

    double dfSdciS[2];
    // 1st species
    dfSdciS[0] = fA_[1] - fC_[1] + derivFMixS[0] + derivIdealMixS[0];

    // 2nd species
    dfSdciS[1] = fB_[1] - fC_[1] + derivFMixS[1] + derivIdealMixS[1];

    // equation fL-fS-(cL-cS)*dfL/dcL=0
    const double fL = cL[0] * fA_[0] + cL[1] * fB_[0]
                      + (1.0 - cL[0] - cL[1]) * fC_[0]
                      + CALPHADcomputeFIdealMixTernary(RT_, cL[0], cL[1])
                      + CALPHADcomputeFMixTernary(
                            L_AB_L_, L_AC_L_, L_BC_L_, L_ABC_L_, cL[0], cL[1]);
    const double fS = cS[0] * fA_[1] + cS[1] * fB_[1]
                      + (1.0 - cS[0] - cS[1]) * fC_[1]
                      + CALPHADcomputeFIdealMixTernary(RT_, cS[0], cS[1])
                      + CALPHADcomputeFMixTernary(
                            L_AB_S_, L_AC_S_, L_BC_S_, L_ABC_S_, cS[0], cS[1]);
    // std::cout<<"fL="<<fL<<", fS="<<fS<<endl;
    fvec[0]
        = fL - fS - (cL[0] - cS[0]) * dfLdciL[0] - (cL[1] - cS[1]) * dfLdciL[1];

    // equation: slope 0 in tangent plane for std::vector orthogonal to cL-CS
    fvec[1] = dfLdciL[0] * cL[1] - dfLdciL[0] * cS[1] - dfLdciL[1] * cL[0]
              + dfLdciL[1] * cS[0];

    fvec[2] = dfLdciL[0] - dfSdciS[0];
    fvec[3] = dfLdciL[1] - dfSdciS[1];

    // std::cout<<"fvec="<<fvec[0]<<","<<fvec[1]<<","<<fvec[2]<<","<<fvec[3]<<endl;
}

//=======================================================================

void CALPHADEqConcentrationSolverTernary::Jacobian(
    const double* const c, double** const fjac)
{
    // tbox::pout<<"Compute Jacobian for CALPHAD..."<<endl;
    const double* const cL = &c[0];
    const double* const cS = &c[2];
    // tbox::pout<<"Compute RHS for CALPHAD..."<<endl;

    double derivIdealMixL[2];
    CALPHADcomputeFIdealMix_derivTernary(RT_, cL[0], cL[1], derivIdealMixL);

    double derivFMixL[2];
    CALPHADcomputeFMix_derivTernary(
        L_AB_L_, L_AC_L_, L_BC_L_, L_ABC_L_, cL[0], cL[1], derivFMixL);

    double deriv2IdealMixL[4];
    CALPHADcomputeFIdealMix_deriv2Ternary(RT_, cL[0], cL[1], deriv2IdealMixL);

    double deriv2FMixL[4];
    CALPHADcomputeFMix_deriv2Ternary(
        L_AB_L_, L_AC_L_, L_BC_L_, L_ABC_L_, cL[0], cL[1], deriv2FMixL);

    double dfLdciL[2];
    // 1st species
    dfLdciL[0] = fA_[0] - fC_[0] + derivFMixL[0] + derivIdealMixL[0];

    // 2nd species
    dfLdciL[1] = fB_[0] - fC_[0] + derivFMixL[1] + derivIdealMixL[1];

    double d2fLdciL2[3]; // include only one cross term (other one equal by
                         // symmetry)
    d2fLdciL2[0] = deriv2FMixL[0] + deriv2IdealMixL[0];
    d2fLdciL2[1] = deriv2FMixL[1] + deriv2IdealMixL[1];
    d2fLdciL2[2] = deriv2FMixL[3] + deriv2IdealMixL[3];

    double derivIdealMixS[2];
    CALPHADcomputeFIdealMix_derivTernary(RT_, cS[0], cS[1], derivIdealMixS);

    double derivFMixS[2];
    CALPHADcomputeFMix_derivTernary(
        L_AB_S_, L_AC_S_, L_BC_S_, L_ABC_S_, cS[0], cS[1], derivFMixS);

    double deriv2IdealMixS[4];
    CALPHADcomputeFIdealMix_deriv2Ternary(RT_, cS[0], cS[1], deriv2IdealMixS);

    double deriv2FMixS[4];
    CALPHADcomputeFMix_deriv2Ternary(
        L_AB_S_, L_AC_S_, L_BC_S_, L_ABC_S_, cS[0], cS[1], deriv2FMixS);

    double d2fSdciS2[3];
    d2fSdciS2[0] = deriv2FMixS[0] + deriv2IdealMixS[0];
    d2fSdciS2[1] = deriv2FMixS[1] + deriv2IdealMixS[1];
    d2fSdciS2[2] = deriv2FMixS[3] + deriv2IdealMixS[3];

    // f[i][j]=df[i]/dc[j]
    fjac[0][0] = fA_[0] - fC_[0] + derivIdealMixL[0] + derivFMixL[0]
                 - dfLdciL[0] - (cL[0] - cS[0]) * d2fLdciL2[0]
                 - (cL[1] - cS[1]) * d2fLdciL2[1];

    fjac[0][1] = fB_[0] - fC_[0] + derivIdealMixL[1] + derivFMixL[1]
                 - dfLdciL[1] - (cL[0] - cS[0]) * d2fLdciL2[1]
                 - (cL[1] - cS[1]) * d2fLdciL2[2];

    fjac[0][2]
        = -fA_[1] + fC_[1] - derivIdealMixS[0] - derivFMixS[0] + dfLdciL[0];

    fjac[0][3]
        = -fB_[1] + fC_[1] - derivIdealMixS[1] - derivFMixS[1] + dfLdciL[1];

    fjac[1][0]
        = d2fLdciL2[0] * (cL[1] - cS[1]) - dfLdciL[1] + d2fLdciL2[1] * cS[0];

    fjac[1][1]
        = dfLdciL[0] - d2fLdciL2[1] * cS[1] + d2fLdciL2[1] * (cL[0] - cS[0]);

    fjac[1][2] = dfLdciL[1];
    fjac[1][3] = -dfLdciL[0];

    fjac[2][0] = d2fLdciL2[0];
    fjac[2][1] = d2fLdciL2[1];
    fjac[2][2] = -d2fSdciS2[0];
    fjac[2][3] = -d2fSdciS2[1];

    fjac[3][0] = d2fLdciL2[1];
    fjac[3][1] = d2fLdciL2[2];
    fjac[3][2] = -d2fSdciS2[1];
    fjac[3][3] = -d2fSdciS2[2];

    // std::cout<<"Jacobian:"<<endl;
    // std::cout<<"("<<fjac[0][0]<<","<<fjac[0][1]<<","<<fjac[0][2]<<","<<fjac[0][3]<<")"<<endl;
    // std::cout<<"("<<fjac[1][0]<<","<<fjac[1][1]<<","<<fjac[1][2]<<","<<fjac[1][3]<<")"<<endl;
    // std::cout<<"("<<fjac[2][0]<<","<<fjac[2][1]<<","<<fjac[2][2]<<","<<fjac[2][3]<<")"<<endl;
    // std::cout<<"("<<fjac[3][0]<<","<<fjac[3][1]<<","<<fjac[3][2]<<","<<fjac[3][3]<<")"<<endl;
}

//=======================================================================

int CALPHADEqConcentrationSolverTernary::ComputeConcentration(
    double* const conc, const double RTinv, const double* const L_AB_L,
    const double* const L_AC_L, const double* const L_BC_L,
    const double* const L_AB_S, const double* const L_AC_S,
    const double* const L_BC_S, const double* const L_ABC_L,
    const double* const L_ABC_S, const double* const fA, const double* const fB,
    const double* const fC)
{
    RTinv_ = RTinv;
    RT_    = 1. / RTinv;

    for (int ii = 0; ii < 4; ii++)
    {
        L_AB_L_[ii] = L_AB_L[ii];
        L_AC_L_[ii] = L_AC_L[ii];
        L_BC_L_[ii] = L_BC_L[ii];
    }
    for (int ii = 0; ii < 4; ii++)
    {
        L_AB_S_[ii] = L_AB_S[ii];
        L_AC_S_[ii] = L_AC_S[ii];
        L_BC_S_[ii] = L_BC_S[ii];
    }
    for (int ii = 0; ii < 3; ii++)
    {
        L_ABC_L_[ii] = L_ABC_L[ii];
        L_ABC_S_[ii] = L_ABC_S[ii];
    }

    // loop over phases (L and S)
    for (int ii = 0; ii < 2; ii++)
    {
        fA_[ii] = fA[ii];
        fB_[ii] = fB[ii];
        fC_[ii] = fC[ii];
    }

    return DampedNewtonSolver::ComputeSolution(conc);
}
}
