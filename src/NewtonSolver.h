#ifndef included_NewtonSolver
#define included_NewtonSolver

namespace Thermo4PFM
{
template <unsigned int Dimension, typename SolverType>
class NewtonSolver
{
public:
    NewtonSolver();

    ~NewtonSolver(){};

    int ComputeSolution(double* conc);

    void SetTolerance(const double t) { tolerance_ = t; }

    void SetMaxIterations(const int m) { max_iters_ = m; }

    void SetVerbose(const bool verbose) { verbose_ = verbose; }

    void SetDamping(const double alpha) { alpha_ = alpha; }

private:
    void UpdateSolution(
        double* const x, const double* const fvec, double** const fjac);

    void internalRHS(const double* const x, double* const fvec)
    {
        static_cast<SolverType*>(this)->RHS(x, fvec);
    }

    double Determinant(double** const mat);
    void CopyMatrix(double** const dst, double** const src);

    void internalJacobian(const double* const x, double** const fjac)
    {
        static_cast<SolverType*>(this)->Jacobian(x, fjac);
    }

    bool CheckTolerance(const double* const fvec);

    int max_iters_;

    // damping factor
    double alpha_;

    double tolerance_;
    bool verbose_;

    // function to compute determinant of matrix of size Dimension
    double (*det_fun_ptr_)(double** const matrix);
};
}

#endif
