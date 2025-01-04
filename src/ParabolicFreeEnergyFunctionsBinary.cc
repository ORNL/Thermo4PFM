#include "ParabolicFreeEnergyFunctionsBinary.h"
#include "ParabolicConcSolverBinary.h"
#include "ParabolicEqConcSolverBinary.h"
#include "functions.h"

#include <iomanip>
#include <string>

namespace Thermo4PFM
{

ParabolicFreeEnergyFunctionsBinary::ParabolicFreeEnergyFunctionsBinary(
    const double Tref, const double coeffL[][2], const double coeffA[][2],
    const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type)
    : Tref_(Tref),
      energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type)
{
    coeffL_[0][0] = coeffL[0][0];
    coeffL_[0][1] = coeffL[0][1];
    coeffL_[1][0] = coeffL[1][0];
    coeffL_[1][1] = coeffL[1][1];
    coeffL_[2][0] = coeffL[2][0];
    coeffL_[2][1] = coeffL[2][1];

    coeffA_[0][0] = coeffA[0][0];
    coeffA_[0][1] = coeffA[0][1];
    coeffA_[1][0] = coeffA[1][0];
    coeffA_[1][1] = coeffA[1][1];
    coeffA_[2][0] = coeffA[2][0];
    coeffA_[2][1] = coeffA[2][1];

    std::string fenergy_diag_filename("energy.vtk");
    fenergy_diag_filename_ = new char[fenergy_diag_filename.length() + 1];
    strcpy(fenergy_diag_filename_, fenergy_diag_filename.c_str());
}

//-----------------------------------------------------------------------
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
void ParabolicFreeEnergyFunctionsBinary::getPhaseCoeffs(const PhaseIndex pi,
    double& a0, double& a1, double& b0, double& b1, double& c0, double& c1)
{
    switch (pi)
    {
        case PhaseIndex::phaseL:
            a0 = coeffL_[0][0];
            a1 = coeffL_[0][1];
            b0 = coeffL_[1][0];
            b1 = coeffL_[1][1];
            c0 = coeffL_[2][0];
            c1 = coeffL_[2][1];
            break;
        case PhaseIndex::phaseA:
            a0 = coeffA_[0][0];
            a1 = coeffA_[0][1];
            b0 = coeffA_[1][0];
            b1 = coeffA_[1][1];
            c0 = coeffA_[2][0];
            c1 = coeffA_[2][1];
            break;
        default:
            break;
    }
}

double ParabolicFreeEnergyFunctionsBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    double a0 = NAN;
    double a1 = NAN;
    double b0 = NAN;
    double b1 = NAN;
    double c0 = NAN;
    double c1 = NAN;
    getPhaseCoeffs(pi, a0, a1, b0, b1, c0, c1);

    const double c = conc[0];
    const double t = (temperature - Tref_);
    double fe = 0.5 * (a1 * t + a0) * c * c + (b1 * t + b0) * c + (c1 * t + c0);

    // subtract -mu*c to get grand potential
    if (gp)
    {
        double deriv;
        computeDerivFreeEnergy(temperature, conc, pi, &deriv);
        fe -= deriv * conc[0];
    }
    return fe;
}

//=======================================================================

void ParabolicFreeEnergyFunctionsBinary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    double a0, a1, b0, b1;
    switch (pi)
    {
        case PhaseIndex::phaseL:
            a0 = coeffL_[0][0];
            a1 = coeffL_[0][1];
            b0 = coeffL_[1][0];
            b1 = coeffL_[1][1];
            break;
        case PhaseIndex::phaseA:
            a0 = coeffA_[0][0];
            a1 = coeffA_[0][1];
            b0 = coeffA_[1][0];
            b1 = coeffA_[1][1];
            break;
        default:
#ifndef HAVE_OPENMP_OFFLOAD
            abort();
#endif
    }

    const double c = conc[0];
    const double t = (temperature - Tref_);
    double mu      = (a1 * t + a0) * c + (b1 * t + b0);

    deriv[0] = mu;
}

//=======================================================================

void ParabolicFreeEnergyFunctionsBinary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    double* d2fdc2)
{
    (void)conc;

    double a0, a1;
    switch (pi)
    {
        case PhaseIndex::phaseL:
            a0 = coeffL_[0][0];
            a1 = coeffL_[0][1];
            break;
        case PhaseIndex::phaseA:
            a0 = coeffA_[0][0];
            a1 = coeffA_[0][1];
            break;
        default:
#ifndef HAVE_OPENMP_OFFLOAD
            abort();
#endif
    }

    d2fdc2[0] = (a1 * (temp - Tref_) + a0);
}

//=======================================================================

bool ParabolicFreeEnergyFunctionsBinary::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
    ParabolicEqConcSolverBinary eq_solver;

    eq_solver.setup(temperature, coeffL_, coeffA_);
    const double newton_tol = 1.e-8;
    int ret = eq_solver.ComputeConcentration(ceq, newton_tol, maxits);

    return (ret >= 0);
}

//=======================================================================

void ParabolicFreeEnergyFunctionsBinary::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc,
    double& fl, double& fa)
{
    // std::cout<<"ParabolicFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;

    double c[2] = { conc, conc };

    ParabolicConcSolverBinary solver;
    solver.setup(conc, hphi[0], hphi[1], temperature - Tref_, coeffL_, coeffA_);
    int ret = solver.ComputeConcentration(c);
    if (ret < 0)
    {
#if 0
        std::cerr << "ERROR in "
                     "ParabolicFreeEnergyFunctionsBinary::"
                     "computePhasesFreeEnergies()"
                     " ---"
                  << "conc=" << conc << ", hphi=" << hphi[0] << std::endl;
        abort();
#endif
    }

    // assert(c[0] >= 0.);
    fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);

    // assert(c[1] >= 0.);
    fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------

int ParabolicFreeEnergyFunctionsBinary::computePhaseConcentrations(
    const double temperature, const double* const conc,
    const double* const hphi, double* x)
{
    // solve system of equations to find (cl,cs) given conc[0] and hphi
    // x: initial guess and solution
    ParabolicConcSolverBinary solver;
    solver.setup(
        conc[0], hphi[0], hphi[1], temperature - Tref_, coeffL_, coeffA_);
    int ret = solver.ComputeConcentration(x);
#if 0
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "ParabolicFreeEnergyFunctionsBinary::"
                     "computePhaseConcentrations() "
                     "failed for conc="
                  << conc[0] << ", hphi0=" << hphi0 << ", hphi1=" << hphi1 << std::endl;
        abort();
    }
#endif

    return ret;
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

//=======================================================================
// compute free energy in [J/mol]

double ParabolicFreeEnergyFunctionsBinary::fchem(
    const double* const phi, const double* const conc, const double temperature)
{
    double hcphi[3];
    hcphi[0] = interp_func(conc_interp_func_type_, phi[0]);
    hcphi[1] = interp_func(conc_interp_func_type_, phi[1]);

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;

    if ((1.0 - phi[0] > tol) && (1.0 - phi[1] > tol))
    {
        computePhasesFreeEnergies(temperature, hcphi, conc[0], fl, fa);
    }
    else
    {
        if (1.0 - phi[0] <= tol)
        {
            fl = computeFreeEnergy(temperature, conc, PhaseIndex::phaseL);
        }
        else
        {
            fa = computeFreeEnergy(temperature, conc, PhaseIndex::phaseA);
        }
    }

    double hfphi[2];
    hfphi[0] = interp_func(energy_interp_func_type_, phi[0]);
    hfphi[1] = interp_func(energy_interp_func_type_, phi[1]);

    return hfphi[0] * fl + hfphi[1] * fa;
}

//=======================================================================

void ParabolicFreeEnergyFunctionsBinary::printEnergyVsComposition(
    const double temperature, std::ostream& os, const double cmin,
    const double cmax, const int npts)
{
    const double dc      = (cmax - cmin) / (double)(npts - 1);
    const double phil[2] = { 1., 0. };
    const double phia[2] = { 0., 1. };

    os << "c, fL, fA" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        double el = fchem(phil, &conc, temperature);
        double ea = fchem(phia, &conc, temperature);

        os << conc << ", " << el << ", " << ea << std::endl;
    }
}
}
