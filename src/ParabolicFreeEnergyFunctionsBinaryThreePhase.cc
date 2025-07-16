#include "ParabolicFreeEnergyFunctionsBinaryThreePhase.h"
#include "ParabolicConcSolverBinaryThreePhase.h"
#include "functions.h"

#include <iomanip>
#include <string>

namespace Thermo4PFM
{

ParabolicFreeEnergyFunctionsBinaryThreePhase::
    ParabolicFreeEnergyFunctionsBinaryThreePhase(const double Tref,
        const double coeffL[][2], const double coeffA[][2],
        const double coeffB[][2],
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type)
    : Tref_(Tref),
      energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type)
{
    aL_[0] = coeffL[0][0];
    aL_[1] = coeffL[0][1];
    bL_[0] = coeffL[1][0];
    bL_[1] = coeffL[1][1];
    cL_[0] = coeffL[2][0];
    cL_[1] = coeffL[2][1];

    aA_[0] = coeffA[0][0];
    aA_[1] = coeffA[0][1];
    bA_[0] = coeffA[1][0];
    bA_[1] = coeffA[1][1];
    cA_[0] = coeffA[2][0];
    cA_[1] = coeffA[2][1];

    aB_[0] = coeffB[0][0];
    aB_[1] = coeffB[0][1];
    bB_[0] = coeffB[1][0];
    bB_[1] = coeffB[1][1];
    cB_[0] = coeffB[2][0];
    cB_[1] = coeffB[2][1];

    std::string fenergy_diag_filename("energy.vtk");
    fenergy_diag_filename_ = new char[fenergy_diag_filename.length() + 1];
    strcpy(fenergy_diag_filename_, fenergy_diag_filename.c_str());
}

//-----------------------------------------------------------------------
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double ParabolicFreeEnergyFunctionsBinaryThreePhase::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    double a0, a1, b0, b1, c0, c1;
    switch (pi)
    {
        case PhaseIndex::phaseL:
            a0 = aL_[0];
            a1 = aL_[1];
            b0 = bL_[0];
            b1 = bL_[1];
            c0 = cL_[0];
            c1 = cL_[1];
            break;
        case PhaseIndex::phaseA:
            a0 = aA_[0];
            a1 = aA_[1];
            b0 = bA_[0];
            b1 = bA_[1];
            c0 = cA_[0];
            c1 = cA_[1];
            break;
        case PhaseIndex::phaseB:
            a0 = aB_[0];
            a1 = aB_[1];
            b0 = bB_[0];
            b1 = bB_[1];
            c0 = cB_[0];
            c1 = cB_[1];
            break;
        default:
            return 0.;
    }

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

void ParabolicFreeEnergyFunctionsBinaryThreePhase::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    double a0, a1, b0, b1;
    switch (pi)
    {
        case PhaseIndex::phaseL:
            a0 = aL_[0];
            a1 = aL_[1];
            b0 = bL_[0];
            b1 = bL_[1];
            break;
        case PhaseIndex::phaseA:
            a0 = aA_[0];
            a1 = aA_[1];
            b0 = bA_[0];
            b1 = bA_[1];
            break;
        case PhaseIndex::phaseB:
            a0 = aB_[0];
            a1 = aB_[1];
            b0 = bB_[0];
            b1 = bB_[1];
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

void ParabolicFreeEnergyFunctionsBinaryThreePhase::
    computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2)
{
    (void)conc;

    double a0, a1;
    switch (pi)
    {
        case PhaseIndex::phaseL:
            a0 = aL_[0];
            a1 = aL_[1];
            break;
        case PhaseIndex::phaseA:
            a0 = aA_[0];
            a1 = aA_[1];
            break;
        case PhaseIndex::phaseB:
            a0 = aB_[0];
            a1 = aB_[1];
            break;
        default:
#ifndef HAVE_OPENMP_OFFLOAD
            abort();
#endif
    }

    d2fdc2[0] = (a1 * (temp - Tref_) + a0);
}

//=======================================================================

void ParabolicFreeEnergyFunctionsBinaryThreePhase::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc,
    double& fl, double& fa, double& fb)
{
    // std::cout<<"ParabolicFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;

    double c[3] = { conc, conc, conc };

    double coeffL[3][2]
        = { { aL_[0], aL_[1] }, { bL_[0], bL_[1] }, { cL_[0], cL_[1] } };
    double coeffA[3][2]
        = { { aA_[0], aA_[1] }, { bA_[0], bA_[1] }, { cA_[0], cA_[1] } };
    double coeffB[3][2]
        = { { aB_[0], aB_[1] }, { bB_[0], bB_[1] }, { cB_[0], cB_[1] } };
    ParabolicConcSolverBinaryThreePhase solver;
    solver.setup(conc, hphi[0], hphi[1], hphi[2], temperature - Tref_, coeffL,
        coeffA, coeffB);
    int ret = solver.ComputeConcentration(c);
    if (ret < 0)
    {
#if 0
        std::cerr << "ERROR in "
                     "ParabolicFreeEnergyFunctionsBinaryThreePhase::"
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

    // assert(c[1] >= 0.);
    fb = computeFreeEnergy(temperature, &c[2], PhaseIndex::phaseB, false);
}

//-----------------------------------------------------------------------

int ParabolicFreeEnergyFunctionsBinaryThreePhase::computePhaseConcentrations(
    const double temperature, const double* const conc, const double* const phi,
    double* x)
{
    const double hphi0 = interp_func(conc_interp_func_type_, phi[0]);
    const double hphi1 = interp_func(conc_interp_func_type_, phi[1]);
    const double hphi2 = interp_func(conc_interp_func_type_, phi[2]);

    // solve system of equations to find (cl,cs) given conc[0] and hphi
    // x: initial guess and solution
    double coeffL[3][2]
        = { { aL_[0], aL_[1] }, { bL_[0], bL_[1] }, { cL_[0], cL_[1] } };
    double coeffA[3][2]
        = { { aA_[0], aA_[1] }, { bA_[0], bA_[1] }, { cA_[0], cA_[1] } };
    double coeffB[3][2]
        = { { aB_[0], aB_[1] }, { bB_[0], bB_[1] }, { cB_[0], cB_[1] } };

    ParabolicConcSolverBinaryThreePhase solver;
    solver.setup(conc[0], hphi0, hphi1, hphi2, temperature - Tref_, coeffL,
        coeffA, coeffB);
    int ret = solver.ComputeConcentration(x);
#if 0
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "ParabolicFreeEnergyFunctionsBinary::"
                     "computePhaseConcentrations() "
                     "failed for conc="
                  << conc[0] << ", hphi0=" << hphi0 << ", hphi1=" << hphi1 << ", hphi2=" << hphi2 << std::endl;
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

double ParabolicFreeEnergyFunctionsBinaryThreePhase::fchem(
    const double* const phi, const double* const conc, const double temperature)
{
    double hcphi[3];
    hcphi[0] = interp_func(conc_interp_func_type_, phi[0]);
    hcphi[1] = interp_func(conc_interp_func_type_, phi[1]);
    hcphi[2] = interp_func(conc_interp_func_type_, phi[2]);

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;
    double fb        = 0.;

    if ((1.0 - phi[0] > tol) && (1.0 - phi[1] > tol) && (1.0 - phi[2] > tol))
    {
        computePhasesFreeEnergies(temperature, hcphi, conc[0], fl, fa, fb);
    }
    else
    {
        if (1.0 - phi[0] <= tol)
        {
            fl = computeFreeEnergy(temperature, conc, PhaseIndex::phaseL);
        }
        else if (1.0 - phi[1] <= tol)
        {
            fa = computeFreeEnergy(temperature, conc, PhaseIndex::phaseA);
        }
        else
        {
            fb = computeFreeEnergy(temperature, conc, PhaseIndex::phaseB);
        }
    }

    double hfphi[3];
    hfphi[0] = interp_func(energy_interp_func_type_, phi[0]);
    hfphi[1] = interp_func(energy_interp_func_type_, phi[1]);
    hfphi[2] = interp_func(energy_interp_func_type_, phi[2]);

    return hfphi[0] * fl + hfphi[1] * fa + hfphi[2] * fb;
}

//=======================================================================

void ParabolicFreeEnergyFunctionsBinaryThreePhase::printEnergyVsComposition(
    const double temperature, std::ostream& os, const double cmin,
    const double cmax, const int npts)
{
    const double dc      = (cmax - cmin) / (double)(npts - 1);
    const double phil[3] = { 1., 0., 0. };
    const double phia[3] = { 0., 1., 0. };
    const double phib[3] = { 0., 0., 1. };

    os << "c, fL, fA, fB" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        double el = fchem(phil, &conc, temperature);
        double ea = fchem(phia, &conc, temperature);
        double eb = fchem(phib, &conc, temperature);

        os << conc << ", " << el << ", " << ea << ", " << eb << std::endl;
    }
}
}
