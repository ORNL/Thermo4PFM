#ifndef Thermo4PFM_included_LinearSolver
#define Thermo4PFM_included_LinearSolver

namespace Thermo4PFM
{
template <unsigned int Dimension, class SolverType, typename JacobianDataType>
class LinearSolver
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// Solve system of equations for tolerance tol, using at most
    /// max_iters iterations
    /// Solution: conc
    int ComputeSolution(double* conc) { return ComputeSolutionInternal(conc); }
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

private:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    void UpdateSolution(double* const x, const double* const fvec,
        JacobianDataType** const fjac);

    int ComputeSolutionInternal(double* const conc);

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

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
};
}

#endif
