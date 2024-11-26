#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"
#include "CALPHADConcSolverBinary.h"
#include "CALPHADEqConcSolverBinary.h"
#include "PhysicalConstants.h"

namespace Thermo4PFM
{

CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB::
    CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB(const double cB,
        boost::property_tree::ptree& input_db,
        boost::optional<boost::property_tree::ptree&> newton_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type)
    : CALPHADFreeEnergyFunctionsBinaryThreePhase(
          input_db, newton_db, energy_interp_func_type, conc_interp_func_type),
      cB_(cB)
{
}

int CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB::
    computePhaseConcentrations(const double temperature, const double* conc,
        const double* const phi, double* x)
{
    const double RT = GASCONSTANT_R_JPKPMOL * temperature;

    CalphadDataType fA[3];
    CalphadDataType fB[3];
    CalphadDataType Lmix_L[4];
    CalphadDataType Lmix_A[4];
    CalphadDataType Lmix_B[4];

    computeTdependentParameters(temperature, Lmix_L, Lmix_A, Lmix_B, fA, fB);

    const double hphi0 = interp_func(conc_interp_func_type_, phi[0]);
    const double hphi1 = interp_func(conc_interp_func_type_, phi[1]);
    const double hphi2 = interp_func(conc_interp_func_type_, phi[2]);

    const double c0 = conc[0] - hphi2 * cB_;

    // solve system of equations to find (cl,cs) given conc[0] and hphi
    // x: initial guess and solution
    CALPHADConcSolverBinary solver;
    // phi = phiA = hphi1
    // 1-phi = phiL = hphi0
    solver.setup(c0, hphi1, hphi0, RT, Lmix_L, Lmix_A, fA, fB);
    int ret = solver.ComputeConcentration(
        x, newton_tol_, newton_maxits_, newton_alpha_);

    x[2] = cB_;

    return ret;
}

bool CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
    // evaluate temperature dependent parameters
    CalphadDataType fA[3];
    CalphadDataType fB[3];

    CalphadDataType Lmix_L[4];
    CalphadDataType Lmix_A[4];
    CalphadDataType Lmix_B[4];

    computeTdependentParameters(temperature, Lmix_L, Lmix_A, Lmix_B, fA, fB);

    double RT = GASCONSTANT_R_JPKPMOL * temperature;

    CALPHADEqConcSolverBinary eq_solver;
    eq_solver.setup(RT, Lmix_L, Lmix_A, fA, fB);
    int ret = eq_solver.ComputeConcentration(ceq, newton_tol_, maxits);

    ceq[2] = cB_;

    return (ret >= 0);
}
}
