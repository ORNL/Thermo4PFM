#include "NewtonSolver.h"
#include "math_utilities.h"

#include <cassert>
#include <cmath>

#define DEBUG_CONVERGENCE
#ifdef DEBUG_CONVERGENCE
#include <iomanip>
#include <iostream>
#include <vector>
#endif

namespace Thermo4PFM
{

NewtonSolver::NewtonSolver(const int ndim)
    : ndim_(ndim),
      max_iters_(50),
      alpha_(1.),
      tolerance_(1.0e-8),
      verbose_(false)
{
    switch (ndim_)
    {
        case 5:
            det_fun_ptr_ = Determinant5;
            break;
        case 4:
            det_fun_ptr_ = Determinant4;
            break;
        case 3:
            det_fun_ptr_ = Determinant3;
            break;
        case 2:
            det_fun_ptr_ = Determinant2;
            break;
        default:
            det_fun_ptr_ = nullptr;
    }
};

//=======================================================================

bool NewtonSolver::CheckTolerance(const double* const fvec)
{
    for (int ii = 0; ii < ndim_; ii++)
    {
        if (fabs(fvec[ii]) >= tolerance_) return false;
    }
    return true;
}

//=======================================================================

void NewtonSolver::CopyMatrix(double** const dst, double** const src)
{
    assert(src != nullptr);
    assert(dst != nullptr);

    for (int jj = 0; jj < ndim_; jj++)
    {
        for (int ii = 0; ii < ndim_; ii++)
        {
            dst[jj][ii] = src[jj][ii];
        }
    }
}

//=======================================================================

double NewtonSolver::Determinant(double** const matrix)
{
    assert(ndim_ == 2 || ndim_ == 3 || ndim_ == 4 || ndim_ == 5);

    return (*det_fun_ptr_)(matrix);
}

//=======================================================================

void NewtonSolver::UpdateSolution(
    double* const c, const double* const fvec, double** const fjac)
{
    double* mem   = new double[ndim_ * (1 + ndim_)];
    double* mtmp  = mem;
    double* del_c = mem + ndim_ * ndim_;

    double** mwork = new double*[ndim_];
    for (int ii = 0; ii < ndim_; ii++)
    {
        mwork[ii] = &mtmp[ii * ndim_];
    }

    const double D     = Determinant(fjac);
    const double D_inv = 1.0 / D;

    // std::cout << "D = " << D << std::endl;

    // use Cramer's rule to solve linear system
    for (int jj = 0; jj < ndim_; jj++)
    {
        CopyMatrix(mwork, fjac);
        for (int ii = 0; ii < ndim_; ii++)
        {
            mwork[ii][jj] = fvec[ii];
        }

        del_c[jj] = D_inv * Determinant(mwork);

        const double maxdel = 0.25;
        if (fabs(del_c[jj]) > maxdel)
            del_c[jj] = del_c[jj] > 0 ? maxdel : -maxdel;

        // std::cout << "del_c[" << jj << "] = " << del_c[jj] << std::endl;
    }

    for (int ii = 0; ii < ndim_; ii++)
    {
        c[ii] = c[ii] - alpha_ * del_c[ii];
    }

    delete[] mwork;
    delete[] mem;
}

//=======================================================================
// conc: initial guess and output solution
//
// Returns number of iterations used, or -1 if not converged
int NewtonSolver::ComputeSolution(double* const conc)
{
    assert(max_iters_ > 1);

    for (int ii = 0; ii < ndim_; ii++)
        assert(conc[ii] == conc[ii]);

#ifdef DEBUG_CONVERGENCE
    std::vector<double> ctmp;
    ctmp.reserve(40);
    // std::cout<<"NewtonSolver::ComputeSolution(), Initial conc=";
    // for(short i=0;i<N;i++)cout<<conc[i]<<",";
    // std::cout<<endl;
#endif

    double* mem  = new double[ndim_ * (1 + ndim_)];
    double* fvec = mem;
    double* ftmp = mem + ndim_;

    double** fjac = new double*[ndim_];
    for (int ii = 0; ii < ndim_; ii++)
    {
        fjac[ii] = &ftmp[ii * ndim_];
    }

    int iterations = 0;
    bool converged = false;

    while (1)
    {

#ifdef DEBUG_CONVERGENCE
        // for ( int ii = 0; ii < ndim_ ; ii++ )cout<<conc[ii]<<endl;
        // std::cout<<endl;

        for (int ii = 0; ii < ndim_; ii++)
            assert(conc[ii] == conc[ii]);
        for (int ii = 0; ii < ndim_; ii++)
            ctmp.push_back(conc[ii]);
#endif
        RHS(conc, fvec);
#ifdef DEBUG_CONVERGENCE
        for (int ii = 0; ii < ndim_; ii++)
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

#ifdef DEBUG_CONVERGENCE
    if (!converged)
    {
        std::cout << std::setprecision(12);
        std::cout << "Concentration history..." << std::endl;
        for (unsigned j = 0; j < ctmp.size(); j = j + ndim_)
        {
            std::cout << "  conc= ";
            for (int ii = 0; ii < ndim_; ii++)
            {
                std::cout << ctmp[j + ii] << "   ";
            }
            std::cout << std::endl;
        }
        for (int ii = 0; ii < ndim_; ii++)
        {
            std::cout << "  conc[" << ii << "] = " << conc[ii] << std::endl;
        }
        for (int ii = 0; ii < ndim_; ii++)
        {
            std::cout << "  rhs[" << ii << "] = " << fvec[ii] << std::endl;
        }
        std::cerr << iterations << " iterations..." << std::endl;
        std::cerr << "Error: too many iterations in NewtonSolver" << std::endl;
    }
#endif

    delete[] fjac;
    delete[] mem;

    if (!converged) return -1;

    return iterations;
}
}
