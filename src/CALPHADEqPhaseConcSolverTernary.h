#ifndef included_CALPHADEqPhaseConcSolverTernary
#define included_CALPHADEqPhaseConcSolverTernary

#include "NewtonSolver.h"
#include "datatypes.h"

namespace Thermo4PFM
{

class CALPHADEqPhaseConcSolverTernary
    : public NewtonSolver<5, CALPHADEqPhaseConcSolverTernary>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    int ComputeConcentration(double* const conc, const double tol,
        const int max_iters, const double alpha = 1.)
    {
        return NewtonSolver::ComputeSolution(conc, tol, max_iters, alpha);
    }

    void setup(const double c0, const double c1, const double RTinv,
        const double* const L_AB_L, const double* const L_AC_L,
        const double* const L_BC_L, const double* const L_AB_S,
        const double* const L_AC_S, const double* const L_BC_S,
        const double* const L_ABC_L, const double* const L_ABC_S,
        const double* const fA, const double* const fB, const double* const fC);

    void RHS(const double* const x, double* const fvec);

    void Jacobian(const double* const x, double** const fjac);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

private:
    double RTinv_;
    double RT_;

    // energies of 3 species, in two phase each
    CalphadDataType fA_[2];
    CalphadDataType fB_[2];
    CalphadDataType fC_[2];

    // L coefficients for 2 possible phases (L and S)
    CalphadDataType L_AB_L_[4];
    CalphadDataType L_AC_L_[4];
    CalphadDataType L_BC_L_[4];
    CalphadDataType L_ABC_L_[3];

    CalphadDataType L_AB_S_[4];
    CalphadDataType L_AC_S_[4];
    CalphadDataType L_BC_S_[4];
    CalphadDataType L_ABC_S_[3];

    // nominal concentration
    double conc_[2];
};
}
#endif
