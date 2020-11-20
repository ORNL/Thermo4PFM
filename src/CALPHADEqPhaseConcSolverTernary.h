#ifndef included_CALPHADEqPhaseConcSolverTernary
#define included_CALPHADEqPhaseConcSolverTernary

#include "DampedNewtonSolver.h"

namespace Thermo4PFM
{

class CALPHADEqPhaseConcentrationSolverTernary : public DampedNewtonSolver
{
public:
    CALPHADEqPhaseConcentrationSolverTernary(const double c0, const double c1);

    virtual ~CALPHADEqPhaseConcentrationSolverTernary(){};

    int ComputeConcentration(double* const conc, const double RTinv,
        const double* const L_AB_L, const double* const L_AC_L,
        const double* const L_BC_L, const double* const L_AB_S,
        const double* const L_AC_S, const double* const L_BC_S,
        const double* const L_ABC_L, const double* const L_ABC_S,
        const double* const fA, const double* const fB, const double* const fC);

    int ComputeConcentration(double* const conc);

    void setup(const double RTinv, const double* const L_AB_L,
        const double* const L_AC_L, const double* const L_BC_L,
        const double* const L_AB_S, const double* const L_AC_S,
        const double* const L_BC_S, const double* const L_ABC_L,
        const double* const L_ABC_S, const double* const fA,
        const double* const fB, const double* const fC);

    virtual void RHS(const double* const x, double* const fvec);

    virtual void Jacobian(const double* const x, double** const fjac);

private:
    double RTinv_;
    double RT_;
    double hphi_;

    // energies of 3 species, in two phase each
    double fA_[2];
    double fB_[2];
    double fC_[2];

    // L coefficients for 2 possible phases (L and S)
    double L_AB_L_[4];
    double L_AC_L_[4];
    double L_BC_L_[4];
    double L_ABC_L_[3];

    double L_AB_S_[4];
    double L_AC_S_[4];
    double L_BC_S_[4];
    double L_ABC_S_[3];

    // nominal concentration
    double conc_[2];
};
}
#endif
