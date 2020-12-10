#ifndef included_CALPHADConcSolverBinary
#define included_CALPHADConcSolverBinary

#include "DampedNewtonSolver.h"

namespace Thermo4PFM
{

class CALPHADConcSolverBinary : public DampedNewtonSolver
{
public:
    CALPHADConcSolverBinary();

    virtual ~CALPHADConcSolverBinary(){};

    int ComputeConcentration(double* const conc, const double c0,
        const double hphi, const double RTinv, const double* const Lmix_L_,
        const double* const Lmix_A_, const double* const fA,
        const double* const fB);

protected:
    double fA_[2];
    double fB_[2];

    // 4 L coefficients for 2 possible phases (L, A)
    double Lmix_L_[4];
    double Lmix_A_[4];
    double RTinv_;

private:
    void computeXi(const double* const c, double xi[2]) const;

    void computeDxiDc(const double* const c, double dxidc[2]) const;

    // virtual functions inherited from DampedNewtonSolver
    void RHS(const double* const x, double* const fvec);

    void Jacobian(const double* const x, double** const fjac);

    double c0_;
    double hphi_;
};
}
#endif
