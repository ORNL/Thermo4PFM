#include "DampedNewtonSolver.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include <iomanip>

namespace Thermo4PFM
{

DampedNewtonSolver::DampedNewtonSolver() : NewtonSolver(), alpha_(1.){};

//=======================================================================
// note: sizes to accomodate up to ternary alloys
//
// c: solution to be updated
void DampedNewtonSolver::UpdateSolution(
    double* const c, const double* const fvec, double** const fjac)
{
    int nn = size();

    static double* mwork[5];
    static double mtmp[25];
    for (int ii = 0; ii < nn; ii++)
    {
        mwork[ii] = &mtmp[ii * nn];
    }

    const double D = Determinant(fjac);
    assert(fabs(D) > 1.e-15);

    const double D_inv = 1.0 / D;

    // std::cout<<setprecision(12);
    // std::cout << "DampedNewtonSolver::UpdateSolution(), N = "<<nn<<", D = "
    // << D << std::endl;

    static double del[5];

    // use Cramer's rule to solve linear system
    for (int jj = 0; jj < nn; jj++)
    {
        CopyMatrix(mwork, fjac);

        // replace jth column with rhs
        for (int ii = 0; ii < nn; ii++)
        {
            mwork[ii][jj] = fvec[ii];
        }

        const double Dmwork = Determinant(mwork);
        // std::cout << "nn="<<nn<<", Dmwork="<<Dmwork <<endl;
        del[jj] = D_inv * Dmwork;

        const double maxdel = 0.25;
        if (fabs(del[jj]) > maxdel) del[jj] = del[jj] > 0 ? maxdel : -maxdel;

        // std::cout << "del[" << jj << "] = " << del[jj] << std::endl;
    }

    for (int ii = 0; ii < nn; ii++)
    {
        c[ii] = c[ii] - alpha_ * del[ii];
    }
}
}
