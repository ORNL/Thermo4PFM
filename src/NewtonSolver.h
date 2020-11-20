#ifndef included_NewtonSolver
#define included_NewtonSolver

class NewtonSolver
{
public:
    NewtonSolver();

    virtual ~NewtonSolver(){};

    virtual void initialize(){};

    virtual int ComputeSolution(double* const conc, const int N);

    void SetTolerance(const double t) { tolerance_ = t; }

    void SetMaxIterations(const int m) { max_iters_ = m; }

    void SetVerbose(const bool verbose) { verbose_ = verbose; }

    double Determinant(double** const mat);

    void CopyMatrix(double** const dst, double** const src);

    int size() const { return s_N; };

    virtual void UpdateSolution(
        double* const x, const double* const fvec, double** const fjac);

    virtual void RHS(const double* const x, double* const fvec) = 0;

private:
    virtual void Jacobian(const double* const x, double** const fjac) = 0;

    bool CheckTolerance(const double* const fvec);
    bool CheckToleranceFirstEq(const double* const fvec);

    /*
     * Number of equations in system
     */
    static int s_N;

    int max_iters_;
    double tolerance_;
    bool verbose_;
};

#endif
