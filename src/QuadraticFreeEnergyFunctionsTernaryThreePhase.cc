#include "QuadraticFreeEnergyFunctionsTernaryThreePhase.h"
#include "QuadraticConcSolverBinaryThreePhase.h"
#include "datatypes.h"
#include "functions.h"
#include "well_functions.h"

#include <cmath>
#include <iomanip>
#include <string>

namespace Thermo4PFM
{
QuadraticFreeEnergyFunctionsTernaryThreePhase::
    QuadraticFreeEnergyFunctionsTernaryThreePhase(const double Al[2],
        const double ceql[2], const double Aa[2], const double ceqa[2],
        const double Ab[2], const double ceqb[2],
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type)
    : Al_{ Al[0], Al[1] },
      Aa_{ Aa[0], Aa[1] },
      Ab_{ Ab[0], Ab[1] },
      ceql_{ ceql[0], ceql[1] },
      ceqa_{ ceqa[0], ceqa[1] },
      ceqb_{ ceqb[0], ceqb[1] },
      energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type)
{
}

//=======================================================================

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

double QuadraticFreeEnergyFunctionsTernaryThreePhase::computeFreeEnergy(
    const double temperature, const double* conc, const PhaseIndex pi,
    const bool gp)
{
    const double* A;
    const double* ceq;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            A   = Al_;
            ceq = ceql_;
            break;
        case PhaseIndex::phaseA:
            A   = Aa_;
            ceq = ceqa_;
            break;
        case PhaseIndex::phaseB:
            A   = Ab_;
            ceq = ceqb_;
            break;
        default:
            //            std::cout <<
            //            "QuadraticFreeEnergyFunctionsTernaryThreePhase::"
            //                         "computeFreeEnergy(), undefined phase!!!"
            //                      << std::endl;
            // abort();
            return 0.;
    }

    double fe = A[0] * (conc[0] - ceq[0]) * (conc[0] - ceq[0])
                + A[1] * (conc[1] - ceq[1]) * (conc[1] - ceq[1]);

    // subtract -mu*c to get grand potential
    if (gp)
    {
        double deriv[2];
        computeDerivFreeEnergy(temperature, conc, pi, deriv);
        fe -= deriv[0] * conc[0];
        fe -= deriv[1] * conc[1];
    }

    return fe;
}

//=======================================================================

void QuadraticFreeEnergyFunctionsTernaryThreePhase::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    (void)temperature;

    const double* A;
    const double* ceq;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            A   = Al_;
            ceq = ceql_;
            break;
        case PhaseIndex::phaseA:
            A   = Aa_;
            ceq = ceqa_;
            break;
        case PhaseIndex::phaseB:
            A   = Ab_;
            ceq = ceqb_;
            break;
        default:
            return;
    }

    deriv[0] = 2. * A[0] * (conc[0] - ceq[0]);
    deriv[1] = 2. * A[1] * (conc[1] - ceq[1]);
}

//=======================================================================

void QuadraticFreeEnergyFunctionsTernaryThreePhase::
    computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2)
{
    (void)temp;
    (void)conc;

    const double* A;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            A = Al_;
            break;
        case PhaseIndex::phaseA:
            A = Aa_;
            break;
        case PhaseIndex::phaseB:
            A = Ab_;
            break;
        default:
            return;
    }

    d2fdc2[0] = A[0];
    d2fdc2[1] = 0.;
    d2fdc2[2] = 0.;
    d2fdc2[3] = A[1];
}

//=======================================================================

void QuadraticFreeEnergyFunctionsTernaryThreePhase::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc0,
    const double conc1, double& fl, double& fa, double& fb)
{
    // std::cout<<"QuadraticFreeEnergyFunctionsTernaryThreePhase::computePhasesFreeEnergies()"<<endl;

    double cauxilliary0[3] = { conc0, conc0, conc0 };
    double cauxilliary1[3] = { conc1, conc1, conc1 };

    QuadraticConcSolverBinaryThreePhase solver0;
    solver0.setup(conc0, hphi[0], hphi[1], hphi[2], Al_[0], ceql_[0], Aa_[0],
        ceqa_[0], Ab_[0], ceqb_[0]);
    int ret = solver0.ComputeConcentration(cauxilliary0);

    QuadraticConcSolverBinaryThreePhase solver1;
    solver1.setup(conc1, hphi[0], hphi[1], hphi[2], Al_[1], ceql_[1], Aa_[1],
        ceqa_[1], Ab_[1], ceqb_[1]);
    ret = solver1.ComputeConcentration(cauxilliary1);

    double concl[2] = { cauxilliary0[0], cauxilliary1[0] };
    fl = computeFreeEnergy(temperature, &concl[0], PhaseIndex::phaseL, false);

    double conca[2] = { cauxilliary0[1], cauxilliary1[1] };
    fa = computeFreeEnergy(temperature, &conca[0], PhaseIndex::phaseA, false);

    double concb[2] = { cauxilliary0[2], cauxilliary1[2] };
    fb = computeFreeEnergy(temperature, &concb[0], PhaseIndex::phaseB, false);
}

//-----------------------------------------------------------------------
// output: x
int QuadraticFreeEnergyFunctionsTernaryThreePhase::computePhaseConcentrations(
    const double temperature, const double* const conc,
    const double* const hphi, double* x)
{
    (void)temperature;

    // solve two binary problems
    QuadraticConcSolverBinaryThreePhase solver0;
    solver0.setup(conc[0], hphi[0], hphi[1], hphi[2], Al_[0], ceql_[0], Aa_[0],
        ceqa_[0], Ab_[0], ceqb_[0]);
    double cauxilliary0[3];
    int ret = solver0.ComputeConcentration(cauxilliary0);

    QuadraticConcSolverBinaryThreePhase solver1;
    solver1.setup(conc[1], hphi[0], hphi[1], hphi[2], Al_[1], ceql_[1], Aa_[1],
        ceqa_[1], Ab_[1], ceqb_[1]);
    double cauxilliary1[3];
    ret = solver1.ComputeConcentration(cauxilliary1);

    x[0] = cauxilliary0[0];
    x[1] = cauxilliary1[0];
    x[2] = cauxilliary0[1];
    x[3] = cauxilliary1[1];
    x[4] = cauxilliary0[2];
    x[5] = cauxilliary1[2];

    return ret;
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

//=======================================================================
// compute free energy in [J/mol]
double QuadraticFreeEnergyFunctionsTernaryThreePhase::fchem(
    const double* const phi, const double* const conc, const double temperature)
{
    const double conc0 = conc[0];
    const double conc1 = conc[1];

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
        computePhasesFreeEnergies(temperature, hcphi, conc0, conc1, fl, fa, fb);
    }
    else
    {
        // don't solve for phases concentrations, just compute energy
        // in either phase
        double conc[2] = { conc0, conc1 };
        if (1.0 - phi[0] <= tol)
        {
            fl = computeFreeEnergy(temperature, &conc[0], PhaseIndex::phaseL);
        }
        else if (1.0 - phi[1] <= tol)
        {
            fa = computeFreeEnergy(temperature, &conc[0], PhaseIndex::phaseA);
        }
    }

    const double hfphi = interp_func(energy_interp_func_type_, phi[0]);

    return (1.0 - hfphi) * fl + hfphi * fa;
}
}
