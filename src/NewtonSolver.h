#ifndef Thermo4PFM_included_NewtonSolver
#define Thermo4PFM_included_NewtonSolver

namespace Thermo4PFM
{
template <unsigned int Dimension, class SolverType, typename JacobianDataType>
class NewtonSolver
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// Solve system of equations for tolerance tol, using at most
    /// max_iters iterations
    /// Solution: conc
    int ComputeSolution(double* conc, const double tol, const int max_iters,
        const double alpha = 1.)
    {
        return ComputeSolutionInternal(conc, tol, max_iters, alpha);
    }
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

private:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    void UpdateSolution(double* const x, const double* const fvec,
        JacobianDataType** const fjac, const double alpha);

    int ComputeSolutionInternal(double* const conc, const double tol,
        const int max_iters, const double alpha = 1.);

    void internalRHS(const double* const x, double* const fvec)
    {
        static_cast<SolverType*>(this)->RHS(x, fvec);
    }

    template <typename ScalarType>
    void CopyMatrix(ScalarType** const dst, ScalarType** const src);

    void internalJacobian(const double* const x, JacobianDataType** const fjac)
    {
        static_cast<SolverType*>(this)->Jacobian(x, fjac);
    }

    bool CheckTolerance(const double* const fvec, const double tol);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
};
}

#endif
