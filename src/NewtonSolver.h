#ifndef included_NewtonSolver
#define included_NewtonSolver

namespace Thermo4PFM
{
template <unsigned int Dimension, class SolverType>
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
        double** const fjac, const double alpha);

    int ComputeSolutionInternal(double* const conc, const double tol,
        const int max_iters, const double alpha = 1.);

    void internalRHS(const double* const x, double* const fvec)
    {
        static_cast<SolverType*>(this)->RHS(x, fvec);
    }

    void CopyMatrix(double** const dst, double** const src);

    void internalJacobian(const double* const x, double** const fjac)
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
