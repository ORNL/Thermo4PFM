#ifndef included_CALPHADConcSolverBinary
#define included_CALPHADConcSolverBinary

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class CALPHADConcSolverBinary : public NewtonSolver<2, CALPHADConcSolverBinary>
{
public:
    CALPHADConcSolverBinary() : NewtonSolver(){};

    ~CALPHADConcSolverBinary(){};

    // compute "internal" concentrations cL, cS by solving KKK
    // equations
    int ComputeConcentration(double* const conc);

    void setup(const double c0, const double hphi, const double RTinv,
        const double* const Lmix_L_, const double* const Lmix_A_,
        const double* const fA, const double* const fB);

    void RHS(const double* const x, double* const fvec);

    void Jacobian(const double* const x, double** const fjac);

private:
    double c0_;
    double hphi_;
    double RTinv_;

    // 4 L coefficients for 2 possible phases (L, A)
    double Lmix_L_[4];
    double Lmix_A_[4];

    double fA_[2];
    double fB_[2];

    // internal functions to help evaluate RHS and Jacobian
    void computeXi(const double* const c, double xi[2]) const;

    void computeDxiDc(const double* const c, double dxidc[2]) const;
};
}
#endif
