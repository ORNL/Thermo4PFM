#ifndef included_NewtonSolver
#define included_NewtonSolver

namespace Thermo4PFM
{

template <unsigned int Dimension>
class NewtonSolver
{
public:
    NewtonSolver();

    virtual ~NewtonSolver(){};

    virtual int ComputeSolution(double* const conc);

    void SetTolerance(const double t) { tolerance_ = t; }

    void SetMaxIterations(const int m) { max_iters_ = m; }

    void SetVerbose(const bool verbose) { verbose_ = verbose; }

    void SetDamping(const double alpha) { alpha_ = alpha; }

    virtual void UpdateSolution(
        double* const x, const double* const fvec, double** const fjac);

    virtual void RHS(const double* const x, double* const fvec) = 0;

protected:
    double Determinant(double** const mat);
    void CopyMatrix(double** const dst, double** const src);

private:
    virtual void Jacobian(const double* const x, double** const fjac) = 0;

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
