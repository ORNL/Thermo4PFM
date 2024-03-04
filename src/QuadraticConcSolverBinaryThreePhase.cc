#include "QuadraticConcSolverBinaryThreePhase.h"

#include <iostream>

namespace Thermo4PFM
{
//=======================================================================

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

//=======================================================================

// solve for c=(c_L, c_A, c_B)
void QuadraticConcSolverBinaryThreePhase::RHS(
    const double* const c, double* const fvec)
{
    fvec[0] = -c0_ + hphi0_ * c[0] + hphi1_ * c[1] + hphi2_ * c[2];
    fvec[1] = Aa_ * (c[1] - ceqa_) - Al_ * (c[0] - ceql_);
    fvec[2] = Ab_ * (c[2] - ceqb_) - Al_ * (c[0] - ceql_);
}

//=======================================================================

void QuadraticConcSolverBinaryThreePhase::Jacobian(
    const double* const c, JacobianDataType** const fjac)
{
    (void)c;

    fjac[0][0] = hphi0_;
    fjac[0][1] = hphi1_;
    fjac[0][2] = hphi2_;

    fjac[1][0] = -Al_;
    fjac[1][1] = Aa_;
    fjac[1][2] = 0.;

    fjac[2][0] = -Al_;
    fjac[2][1] = 0.;
    fjac[2][2] = Ab_;
}

// set values of internal variables
void QuadraticConcSolverBinaryThreePhase::setup(const double c0,
    const double hphi0, const double hphi1, const double hphi2, const double Al,
    const double ceql, const double Aa, const double ceqa, const double Ab,
    const double ceqb)
{
    c0_    = c0;
    hphi0_ = hphi0;
    hphi1_ = hphi1;
    hphi2_ = hphi2;

    Al_ = Al;
    Aa_ = Aa;
    Ab_ = Ab;

    ceql_ = ceql;
    ceqa_ = ceqa;
    ceqb_ = ceqb;
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
