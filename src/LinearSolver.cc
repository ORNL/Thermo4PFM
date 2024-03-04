#include "QuadraticConcSolverBinaryThreePhase.h"

#include "Determinant.h"
#include "LinearSolver.h"

#ifndef HAVE_OPENMP_OFFLOAD
#include <cassert>
#endif
#include <cmath>
#include <iostream>

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

template <unsigned int Dimension, typename SolverType,
    typename JacobianDataType>
template <typename ScalarType>
void LinearSolver<Dimension, SolverType, JacobianDataType>::CopyMatrix(
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
void LinearSolver<Dimension, SolverType, JacobianDataType>::UpdateSolution(
    double* const c, const double* const fvec, JacobianDataType** const fjac)
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
    }

    for (unsigned int ii = 0; ii < Dimension; ii++)
    {
        c[ii] = c[ii] - del_c[ii];
    }
}

//=======================================================================
// conc: output solution
//
template <unsigned int Dimension, typename SolverType,
    typename JacobianDataType>
int LinearSolver<Dimension, SolverType,
    JacobianDataType>::ComputeSolutionInternal(double* conc)
{
    // for (int ii = 0; ii < Dimension; ii++)
    //    assert(conc[ii] == conc[ii]);

    double fvec[Dimension];

    JacobianDataType ftmp[Dimension * Dimension];
    JacobianDataType* fjac[Dimension];
    for (unsigned int ii = 0; ii < Dimension; ii++)
    {
        fjac[ii] = &ftmp[ii * Dimension];
    }

#ifndef HAVE_OPENMP_OFFLOAD
    for (unsigned int ii = 0; ii < Dimension; ii++)
        assert(conc[ii] == conc[ii]);
#endif
    internalRHS(conc, fvec);
#ifndef HAVE_OPENMP_OFFLOAD
    for (unsigned int ii = 0; ii < Dimension; ii++)
        assert(fvec[ii] == fvec[ii]);
#endif

    internalJacobian(conc, fjac);
    UpdateSolution(conc, fvec, fjac);

    return 0;
}

template class LinearSolver<3, QuadraticConcSolverBinaryThreePhase,
    JacobianDataType>;
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
