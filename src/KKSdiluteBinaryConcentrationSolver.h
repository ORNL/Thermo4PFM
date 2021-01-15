#ifndef included_KKSdiluteBinaryConcentrationSolver
#define included_KKSdiluteBinaryConcentrationSolver

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class KKSdiluteBinaryConcentrationSolver : public NewtonSolver<2>
{
public:
    KKSdiluteBinaryConcentrationSolver();

    ~KKSdiluteBinaryConcentrationSolver(){};

    int ComputeConcentration(double* const conc, const double c0,
        const double hphi, const double RTinv, const double fA,
        const double fB);

private:
    double fA_;
    double fB_;

    void RHS(const double* const x, double* const fvec);

    void Jacobian(const double* const x, double** const fjac);

    double c0_;
    double hphi_;
};
}

#endif
