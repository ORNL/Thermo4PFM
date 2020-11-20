#include "NewtonSolver.h"
#include "math_utilities.h"

#include <cassert>
#include <cmath>
#include <iostream>

#define DEBUG_CONVERGENCE
#ifdef DEBUG_CONVERGENCE
#include <vector>
#endif

#include <iomanip>

namespace Thermo4PFM
{

int NewtonSolver::s_N = 0;

NewtonSolver::NewtonSolver()
    : max_iters_(50), tolerance_(1.0e-8), verbose_(false){};

//=======================================================================

bool NewtonSolver::CheckTolerance(const double* const fvec)
{
    for (int ii = 0; ii < s_N; ii++)
    {
        if (std::abs(fvec[ii]) >= tolerance_) return false;
    }
    return true;
}

//=======================================================================

bool NewtonSolver::CheckToleranceFirstEq(const double* const fvec)
{
    if (std::abs(fvec[0]) >= tolerance_) return false;
    return true;
}

//=======================================================================

void NewtonSolver::CopyMatrix(double** const dst, double** const src)
{
    assert(src != nullptr);
    assert(dst != nullptr);

    for (int jj = 0; jj < s_N; jj++)
    {
        for (int ii = 0; ii < s_N; ii++)
        {
            dst[jj][ii] = src[jj][ii];
        }
    }
}

//=======================================================================

double NewtonSolver::Determinant(double** const m)
{
    assert(s_N == 2 || s_N == 3 || s_N == 4 || s_N == 5);

    if (s_N == 5)
    {
        return DeterminantN(m, 5);
    }
    else if (s_N == 4)
    {
        return Determinant4(m);
    }
    else if (s_N == 3)
    {
        return Determinant3(m);
    }
    else if (s_N == 2)
    {
        return m[0][0] * m[1][1] - m[1][0] * m[0][1];
    }

    return 0.;
}

//=======================================================================

void NewtonSolver::UpdateSolution(
    double* const c, const double* const fvec, double** const fjac)
{
    static double* mwork[3];
    static double mtmp[9];
    for (int ii = 0; ii < s_N; ii++)
    {
        mwork[ii] = &mtmp[ii * s_N];
    }

    const double D     = Determinant(fjac);
    const double D_inv = 1.0 / D;

    // std::cout << "D = " << D << std::endl;

    static double del_c[4];

    // use Cramer's rule to solve linear system
    for (int jj = 0; jj < s_N; jj++)
    {

        CopyMatrix(mwork, fjac);
        for (int ii = 0; ii < s_N; ii++)
        {
            mwork[ii][jj] = fvec[ii];
        }

        del_c[jj] = D_inv * Determinant(mwork);

        // std::cout << "del_c[" << jj << "] = " << del_c[jj] << std::endl;
    }

    double w = 1.0;
    for (int ii = 0; ii < s_N; ii++)
    {
        c[ii] = c[ii] - w * del_c[ii];
    }
}

//=======================================================================
// conc: initial guess and output solution
int NewtonSolver::ComputeSolution(double* const conc, const int N)
{
    assert(max_iters_ > 1);
    for (int ii = 0; ii < N; ii++)
        assert(conc[ii] == conc[ii]);

#ifdef DEBUG_CONVERGENCE
    std::vector<double> ctmp;
    ctmp.reserve(40);
    // std::cout<<"NewtonSolver::ComputeSolution(), Initial conc=";
    // for(short i=0;i<N;i++)cout<<conc[i]<<",";
    // std::cout<<endl;
#endif

    static double* fvec = nullptr;
    static double** fjac;
    static double* ftmp = nullptr;
    if (ftmp == nullptr || s_N != N)
    {
        if (fvec != nullptr) delete[] fvec;
        if (fjac != nullptr) delete[] fjac;
        if (ftmp != nullptr) delete[] ftmp;

        s_N  = N;
        fvec = new double[N];
        fjac = new double*[N];
        ftmp = new double[N * N];
        for (int ii = 0; ii < s_N; ii++)
        {
            fjac[ii] = &ftmp[ii * s_N];
        }
    }

    int iterations = 0;
    bool converged = false;

    initialize();

    while (1)
    {

#ifdef DEBUG_CONVERGENCE
        // for ( int ii = 0; ii < N ; ii++ )cout<<conc[ii]<<endl;
        // std::cout<<endl;

        for (int ii = 0; ii < N; ii++)
            assert(conc[ii] == conc[ii]);
        for (int ii = 0; ii < N; ii++)
            ctmp.push_back(conc[ii]);
#endif
        RHS(conc, fvec);
#ifdef DEBUG_CONVERGENCE
        for (int ii = 0; ii < N; ii++)
            assert(fvec[ii] == fvec[ii]);
#endif

        if (CheckTolerance(fvec))
        {
            converged = true;
            break;
        }

        if (iterations == max_iters_) break;

        Jacobian(conc, fjac);

        UpdateSolution(conc, fvec, fjac);

        iterations++;
    }

    if (!converged)
    {
#ifdef DEBUG_CONVERGENCE
        std::cout << std::setprecision(12);
        std::cout << "Concentration history..." << std::endl;
        for (unsigned j = 0; j < ctmp.size(); j = j + s_N)
        {
            std::cout << "  conc= ";
            for (int ii = 0; ii < s_N; ii++)
            {
                std::cout << ctmp[j + ii] << "   ";
            }
            std::cout << std::endl;
        }
        for (int ii = 0; ii < s_N; ii++)
        {
            std::cout << "  conc[" << ii << "] = " << conc[ii] << std::endl;
        }
        for (int ii = 0; ii < s_N; ii++)
        {
            std::cout << "  rhs[" << ii << "] = " << fvec[ii] << std::endl;
        }
#endif
        std::cerr << iterations << " iterations..." << std::endl;
        std::cerr << "Error: too many iterations in NewtonSolver" << std::endl;
        return -1;
    }

    return iterations;
}
}
