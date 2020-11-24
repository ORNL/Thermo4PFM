#ifndef included_DampedNewtonSolver
#define included_DampedNewtonSolver

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class DampedNewtonSolver : public NewtonSolver
{
public:
    DampedNewtonSolver(const int ndim);

    virtual ~DampedNewtonSolver(){};

    void SetDamping(const double alpha) { alpha_ = alpha; }

    virtual void UpdateSolution(
        double* const x, const double* const fvec, double** const fjac);

private:
    // damping factor
    double alpha_;
};
}
#endif
