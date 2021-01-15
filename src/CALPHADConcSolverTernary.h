#ifndef included_CALPHADConcSolverTernary
#define included_CALPHADConcSolverTernary

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class CALPHADConcentrationSolverTernary : public NewtonSolver<4>
{
public:
    CALPHADConcentrationSolverTernary();

    virtual ~CALPHADConcentrationSolverTernary(){};

    void setup(const double c0, const double c1, const double hphi,
        const double RTinv, const double* const L_AB_L,
        const double* const L_AC_L, const double* const L_BC_L,
        const double* const L_AB_S, const double* const L_AC_S,
        const double* const L_BC_S, const double* const L_ABC_L,
        const double* const L_ABC_S, const double* const fA,
        const double* const fB, const double* const fC);

    int ComputeConcentration(double* const conc);

    // virtual functions inherited from NewtonSolver
    void RHS(const double* const c, double* const fvec);

    void Jacobian(const double* const c, double** const fjac);

private:
    // energies of 3 species, in two phase each
    double fA_[2];
    double fB_[2];
    double fC_[2];

    // L coefficients for phase L
    double L_AB_L_[4];
    double L_AC_L_[4];
    double L_BC_L_[4];
    double L_ABC_L_[3];

    // L coefficients for phase S
    double L_AB_S_[4];
    double L_AC_S_[4];
    double L_BC_S_[4];
    double L_ABC_S_[3];

    double c0_[2];
    double hphi_;

    double RT_;
};
}
#endif
