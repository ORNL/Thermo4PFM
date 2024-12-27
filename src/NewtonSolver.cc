#include "CALPHADConcSolverBinary.h"
#include "CALPHADConcSolverBinary2Ph1Sl.h"
#include "CALPHADConcSolverBinary3Ph2Sl.h"
#include "CALPHADConcSolverBinaryThreePhase.h"
#include "CALPHADConcSolverTernary.h"
#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADEqConcSolverBinary2Ph1Sl.h"
#include "CALPHADEqConcSolverTernary.h"
#include "CALPHADTieLineConcSolverTernary.h"
#include "KKSdiluteBinaryConcSolver.h"
#include "ParabolicEqConcSolverBinary.h"

#include "Determinant.h"
#include "NewtonSolver.h"

#ifndef HAVE_OPENMP_OFFLOAD
#include <cassert>
#endif
#include <cmath>

#ifdef WITH_CONVERGENCE_HISTORY
#include <iomanip>
#include <iostream>
#include <vector>
#endif

#include <iostream>
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
    for (unsigned int ii = 0; ii < Dimension; ii++)
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

    for (unsigned int jj = 0; jj < Dimension; jj++)
    {
        for (unsigned int ii = 0; ii < Dimension; ii++)
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
    for (unsigned int ii = 0; ii < Dimension; ii++)
    {
        mwork[ii] = &memf[ii * Dimension];
    }

    const double D     = evalDeterminant<Dimension, JacobianDataType>(fjac);
    const double D_inv = 1.0 / D;

    // use Cramer's rule to solve linear system
    for (unsigned int jj = 0; jj < Dimension; jj++)
    {
        CopyMatrix(mwork, fjac);
        for (unsigned int ii = 0; ii < Dimension; ii++)
        {
            mwork[ii][jj] = fvec[ii];
        }

        del_c[jj] = D_inv * evalDeterminant<Dimension, JacobianDataType>(mwork);

        const double maxdel = 0.25;
        if (fabs(del_c[jj]) > maxdel)
            del_c[jj] = del_c[jj] > 0 ? maxdel : -maxdel;

        // std::cout << "del_c[" << jj << "] = " << del_c[jj] << std::endl;
    }

    for (unsigned int ii = 0; ii < Dimension; ii++)
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
    // for (unsigned int ii = 0; ii < Dimension; ii++)
    //    assert(conc[ii] == conc[ii]);

#ifdef WITH_CONVERGENCE_HISTORY
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
    for (unsigned int ii = 0; ii < Dimension; ii++)
    {
        fjac[ii] = &ftmp[ii * Dimension];
    }

    int iterations = 0;
    bool converged = false;

    while (1)
    {
#ifndef HAVE_OPENMP_OFFLOAD
        for (unsigned int ii = 0; ii < Dimension; ii++)
            assert(conc[ii] == conc[ii]);
#endif
#ifdef WITH_CONVERGENCE_HISTORY
        for (unsigned int ii = 0; ii < Dimension; ii++)
            ctmp.push_back(conc[ii]);
#endif
        internalRHS(conc, fvec);
#ifndef HAVE_OPENMP_OFFLOAD
        for (unsigned int ii = 0; ii < Dimension; ii++)
            assert(fvec[ii] == fvec[ii]);
#endif
#ifdef WITH_CONVERGENCE_HISTORY
        for (unsigned int ii = 0; ii < Dimension; ii++)
            residual.push_back(fvec[ii]);
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

#ifdef WITH_CONVERGENCE_HISTORY
    std::cout << std::setprecision(12);
    std::cout << "======================" << std::endl;
    std::cout << "Convergence history..." << std::endl;
    for (unsigned j = 0; j < ctmp.size(); j = j + Dimension)
    {
        std::cout << "  conc= ";
        for (unsigned int ii = 0; ii < Dimension; ii++)
        {
            std::cout << ctmp[j + ii] << "   ";
        }
        std::cout << ", rhs = ";
        for (unsigned int ii = 0; ii < Dimension; ii++)
        {
            std::cout << residual[j + ii] << "   ";
        }
        std::cout << std::endl;
    }
    std::cout << "======================" << std::endl;
    std::cout << "Final solution:" << std::endl;
    for (unsigned int ii = 0; ii < Dimension; ii++)
    {
        std::cout << "  conc[" << ii << "] = " << conc[ii] << std::endl;
    }
    for (unsigned int ii = 0; ii < Dimension; ii++)
    {
        std::cout << "  rhs[" << ii << "] = " << fvec[ii] << std::endl;
    }
#endif

    if (!converged)
    {
#ifndef HAVE_OPENMP_OFFLOAD
        // just print a warning message and let caller decide what to do about
        // it
        std::cerr << "WARNING: too many iterations in NewtonSolver"
                  << std::endl;
        std::cerr << iterations << " iterations..." << std::endl;
        for (unsigned int ii = 0; ii < Dimension; ii++)
        {
            std::cout << "  conc[" << ii << "] = " << conc[ii] << "  rhs[" << ii
                      << "] = " << fvec[ii] << std::endl;
        }
#endif
        return -1;
    }
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
template class NewtonSolver<2, CALPHADConcSolverBinary2Ph1Sl, JacobianDataType>;
template class NewtonSolver<2, CALPHADEqConcSolverBinary2Ph1Sl,
    JacobianDataType>;
template class NewtonSolver<2, ParabolicEqConcSolverBinary, JacobianDataType>;
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
