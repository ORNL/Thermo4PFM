#ifndef included_NewtonSolver
#define included_NewtonSolver

namespace Thermo4PFM
{
template <unsigned int Dimension, typename SolverType>
class NewtonSolver
{
public:
    int ComputeSolution(double* conc, const double tol, const int max_iters,
        const double alpha = 1.);

private:
    void UpdateSolution(double* const x, const double* const fvec,
        double** const fjac, const double alpha);

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
};
}

#endif
