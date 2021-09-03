#include "CALPHADConcSolverBinary.h"
#include "CALPHADConcSolverBinary3Ph2Sl.h"
#include "CALPHADConcSolverBinaryThreePhase.h"
#include "CALPHADConcSolverTernary.h"
#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADEqConcSolverTernary.h"
#include "CALPHADTieLineConcSolverTernary.h"
#include "KKSdiluteBinaryConcSolver.h"

#include "Determinant.h"
#include "NewtonSolver.h"

#ifndef HAVE_OPENMP_OFFLOAD
#include <cassert>
#endif
#include <cmath>

//#define DEBUG_CONVERGENCE
#ifdef DEBUG_CONVERGENCE
#include <iomanip>
#include <iostream>
#include <vector>
#endif

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
template <unsigned int Dimension, typename SolverType,
    typename JacobianDataType>
bool NewtonSolver<Dimension, SolverType, JacobianDataType>::CheckTolerance(
    const double* const fvec, const double tol)
{
    bool ret = true;
    for (int ii = 0; ii < Dimension; ii++)
    {
        ret = (ret && (fabs(fvec[ii]) < tol));
    }
    return ret;
}

//=======================================================================

template <unsigned int Dimension, typename SolverType,
    typename JacobianDataType>
template <typename ScalarType>
void NewtonSolver<Dimension, SolverType, JacobianDataType>::CopyMatrix(
    ScalarType** const dst, ScalarType** const src)
{
#ifndef HAVE_OPENMP_OFFLOAD
    assert(src != nullptr);
    assert(dst != nullptr);
#endif

    for (int jj = 0; jj < Dimension; jj++)
    {
        for (int ii = 0; ii < Dimension; ii++)
        {
            dst[jj][ii] = src[jj][ii];
        }
    }
}

//=======================================================================

template <unsigned int Dimension, typename SolverType,
    typename JacobianDataType>
void NewtonSolver<Dimension, SolverType, JacobianDataType>::UpdateSolution(
    double* const c, const double* const fvec, JacobianDataType** const fjac,
    const double alpha)
{
    double del_c[Dimension];

    JacobianDataType memf[Dimension * Dimension];
    JacobianDataType* mwork[Dimension];
    for (int ii = 0; ii < Dimension; ii++)
    {
        mwork[ii] = &memf[ii * Dimension];
    }

    const double D     = evalDeterminant<Dimension, JacobianDataType>(fjac);
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

        del_c[jj] = D_inv * evalDeterminant<Dimension, JacobianDataType>(mwork);

        const double maxdel = 0.25;
        if (fabs(del_c[jj]) > maxdel)
            del_c[jj] = del_c[jj] > 0 ? maxdel : -maxdel;

        // std::cout << "del_c[" << jj << "] = " << del_c[jj] << std::endl;
    }

    for (int ii = 0; ii < Dimension; ii++)
    {
        c[ii] = c[ii] - alpha * del_c[ii];
    }
}

//=======================================================================
// conc: initial guess and output solution
//
// Returns number of iterations used, or -1 if not converged
template <unsigned int Dimension, typename SolverType,
    typename JacobianDataType>
int NewtonSolver<Dimension, SolverType,
    JacobianDataType>::ComputeSolutionInternal(double* conc, const double tol,
    const int max_iters, const double alpha)
{
    // assert(max_iters > 1);
    // for (int ii = 0; ii < Dimension; ii++)
    //    assert(conc[ii] == conc[ii]);

#ifdef DEBUG_CONVERGENCE
    std::vector<double> ctmp;
    ctmp.reserve(40);
    std::vector<double> residual;
    // std::cout<<"NewtonSolver::ComputeSolution(), Initial conc=";
    // for(short i=0;i<N;i++)cout<<conc[i]<<",";
    // std::cout<<endl;
#endif

    double fvec[Dimension];

    JacobianDataType ftmp[Dimension * Dimension];
    JacobianDataType* fjac[Dimension];
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
        internalRHS(conc, fvec);
#ifdef DEBUG_CONVERGENCE
        for (int ii = 0; ii < Dimension; ii++)
            residual.push_back(fvec[ii]);
        for (int ii = 0; ii < Dimension; ii++)
            assert(fvec[ii] == fvec[ii]);
#endif

        if (CheckTolerance(fvec, tol))
        {
            converged = true;
            break;
        }

        if (iterations == max_iters) break;

        internalJacobian(conc, fjac);
        UpdateSolution(conc, fvec, fjac, alpha);

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
            for (int ii = 0; ii < Dimension; ii++)
            {
                std::cout << residual[j + ii] << "   ";
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

template class NewtonSolver<2, CALPHADConcSolverBinary, JacobianDataType>;
template class NewtonSolver<3, CALPHADConcSolverBinaryThreePhase,
    JacobianDataType>;
template class NewtonSolver<3, CALPHADConcSolverBinary3Ph2Sl, JacobianDataType>;
template class NewtonSolver<2, CALPHADEqConcSolverBinary, JacobianDataType>;
template class NewtonSolver<4, CALPHADConcSolverTernary, JacobianDataType>;
template class NewtonSolver<4, CALPHADEqConcSolverTernary, JacobianDataType>;
template class NewtonSolver<5, CALPHADTieLineConcSolverTernary,
    JacobianDataType>;
template class NewtonSolver<2, KKSdiluteBinaryConcSolver, JacobianDataType>;
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
