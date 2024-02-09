#ifndef Thermo4PFM_included_QuadraticConcSolverBinaryThreePhase
#define Thermo4PFM_included_QuadraticConcSolverBinaryThreePhase

#include "LinearSolver.h"
#include "datatypes.h"

namespace Thermo4PFM
{
class QuadraticConcSolverBinaryThreePhase
    : public LinearSolver<3, QuadraticConcSolverBinaryThreePhase,
          JacobianDataType>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// compute "internal" concentrations cL, cS1, cS2 by solving KKS
    /// equations
    /// conc: solution (concentration in each phase)
    int ComputeConcentration(double* const conc)
    {
        return LinearSolver::ComputeSolution(conc);
    }

    /// setup model paramater values to be used by solver,
    /// including composition "c0" and phase fraction "hphi"
    /// to solve for
    void setup(const double c0, const double hphi0, const double hphi1,
        const double hphi2, const double Al, const double ceql, const double Aa,
        const double ceqa, const double Ab, const double ceqb);

    /// evaluate RHS of the system of eqautions to solve for
    /// specific to this solver
    void RHS(const double* const x, double* const fvec);

    /// evaluate Jacobian of system of equations
    /// specific to this solver
    void Jacobian(const double* const x, JacobianDataType** const fjac);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

private:
    ///
    /// Nominal composition to solve for
    ///
    double c0_;

    ///
    /// phase fractions to solve for
    ///
    double hphi0_;
    double hphi1_;
    double hphi2_;

    double scaledRT_;

    ///
    /// factors for quadratic energy for 3 possible phases (L, A, B)
    ///
    double Al_;
    double Aa_;
    double Ab_;

    ///
    /// equilibrium compositions for quadratic energy for 3 possible phases (L,
    /// A,B)
    ///
    double ceql_;
    double ceqa_;
    double ceqb_;
};
}
#endif
