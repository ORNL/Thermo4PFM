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

template <unsigned int Dimension>
NewtonSolver<Dimension>::NewtonSolver()
    : max_iters_(50), alpha_(1.), tolerance_(1.0e-8), verbose_(false)
{
    switch (Dimension)
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

template <unsigned int Dimension>
bool NewtonSolver<Dimension>::CheckTolerance(const double* const fvec)
{
    for (int ii = 0; ii < Dimension; ii++)
    {
        if (fabs(fvec[ii]) >= tolerance_) return false;
    }
    return true;
}

//=======================================================================

template <unsigned int Dimension>
void NewtonSolver<Dimension>::CopyMatrix(double** const dst, double** const src)
{
    assert(src != nullptr);
    assert(dst != nullptr);

    for (int jj = 0; jj < Dimension; jj++)
    {
        for (int ii = 0; ii < Dimension; ii++)
        {
            dst[jj][ii] = src[jj][ii];
        }
    }
}

//=======================================================================

template <unsigned int Dimension>
double NewtonSolver<Dimension>::Determinant(double** const matrix)
{
    return (*det_fun_ptr_)(matrix);
}

//=======================================================================

template <unsigned int Dimension>
void NewtonSolver<Dimension>::UpdateSolution(
    double* const c, const double* const fvec, double** const fjac)
{
    double mem[Dimension * (1 + Dimension)];
    double* mtmp  = mem;
    double* del_c = mem + Dimension * Dimension;

    double* mwork[Dimension];
    for (int ii = 0; ii < Dimension; ii++)
    {
        mwork[ii] = &mtmp[ii * Dimension];
    }

    const double D     = Determinant(fjac);
    const double D_inv = 1.0 / D;

    // std::cout << "D = " << D << std::endl;

    // use Cramer's rule to solve linear system
    for (int jj = 0; jj < Dimension; jj++)
    {
        CopyMatrix(mwork, fjac);
        for (int ii = 0; ii < Dimension; ii++)
        {
            mwork[ii][jj] = fvec[ii];
        }

        del_c[jj] = D_inv * Determinant(mwork);

        const double maxdel = 0.25;
        if (fabs(del_c[jj]) > maxdel)
            del_c[jj] = del_c[jj] > 0 ? maxdel : -maxdel;

        // std::cout << "del_c[" << jj << "] = " << del_c[jj] << std::endl;
    }

    for (int ii = 0; ii < Dimension; ii++)
    {
        c[ii] = c[ii] - alpha_ * del_c[ii];
    }
}

//=======================================================================
// conc: initial guess and output solution
//
// Returns number of iterations used, or -1 if not converged
template <unsigned int Dimension>
int NewtonSolver<Dimension>::ComputeSolution(double* const conc)
{
    assert(max_iters_ > 1);

    for (int ii = 0; ii < Dimension; ii++)
        assert(conc[ii] == conc[ii]);

#ifdef DEBUG_CONVERGENCE
    std::vector<double> ctmp;
    ctmp.reserve(40);
    // std::cout<<"NewtonSolver::ComputeSolution(), Initial conc=";
    // for(short i=0;i<N;i++)cout<<conc[i]<<",";
    // std::cout<<endl;
#endif

    double mem[Dimension * (1 + Dimension)];
    double* fvec = mem;
    double* ftmp = mem + Dimension;

    double* fjac[Dimension];
    for (int ii = 0; ii < Dimension; ii++)
    {
        fjac[ii] = &ftmp[ii * Dimension];
    }

    int iterations = 0;
    bool converged = false;

    while (1)
    {

#ifdef DEBUG_CONVERGENCE
        // for ( int ii = 0; ii < Dimension ; ii++ )cout<<conc[ii]<<endl;
        // std::cout<<endl;

        for (int ii = 0; ii < Dimension; ii++)
            assert(conc[ii] == conc[ii]);
        for (int ii = 0; ii < Dimension; ii++)
            ctmp.push_back(conc[ii]);
#endif
        RHS(conc, fvec);
#ifdef DEBUG_CONVERGENCE
        for (int ii = 0; ii < Dimension; ii++)
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
        for (unsigned j = 0; j < ctmp.size(); j = j + Dimension)
        {
            std::cout << "  conc= ";
            for (int ii = 0; ii < Dimension; ii++)
            {
                std::cout << ctmp[j + ii] << "   ";
            }
            std::cout << std::endl;
        }
        for (int ii = 0; ii < Dimension; ii++)
        {
            std::cout << "  conc[" << ii << "] = " << conc[ii] << std::endl;
        }
        for (int ii = 0; ii < Dimension; ii++)
        {
            std::cout << "  rhs[" << ii << "] = " << fvec[ii] << std::endl;
        }
        std::cerr << iterations << " iterations..." << std::endl;
        std::cerr << "Error: too many iterations in NewtonSolver" << std::endl;
    }
#endif

    if (!converged) return -1;

    return iterations;
}

template class NewtonSolver<2>;
template class NewtonSolver<4>;
template class NewtonSolver<5>;
}
