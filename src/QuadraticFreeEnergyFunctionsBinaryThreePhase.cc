#include "QuadraticFreeEnergyFunctionsBinaryThreePhase.h"
#include "QuadraticConcSolverBinaryThreePhase.h"
#include "functions.h"

#include <iomanip>
#include <string>

namespace Thermo4PFM
{

QuadraticFreeEnergyFunctionsBinaryThreePhase::
    QuadraticFreeEnergyFunctionsBinaryThreePhase(const double Al,
        const double ceql, const double Aa, const double ceqa, const double Ab,
        const double ceqb,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type)
    : Al_(Al),
      ceql_(ceql),
      Aa_(Aa),
      ceqa_(ceqa),
      Ab_(Ab),
      ceqb_(ceqb),
      energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type)
{
    std::string fenergy_diag_filename("energy.vtk");
    fenergy_diag_filename_ = new char[fenergy_diag_filename.length() + 1];
    strcpy(fenergy_diag_filename_, fenergy_diag_filename.c_str());
}

//-----------------------------------------------------------------------
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double QuadraticFreeEnergyFunctionsBinaryThreePhase::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    double A;
    double ceq;

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
            return 0.;
    }

    double fe = A * (conc[0] - ceq) * (conc[0] - ceq);

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

void QuadraticFreeEnergyFunctionsBinaryThreePhase::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    double A;
    double ceq;

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

    double mu = 2. * A * (conc[0] - ceq);

    deriv[0] = mu;
}

//=======================================================================

void QuadraticFreeEnergyFunctionsBinaryThreePhase::
    computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2)
{
    double deriv;
    switch (pi)
    {
        case PhaseIndex::phaseL:
            deriv = 2. * Al_;
            break;
        case PhaseIndex::phaseA:
            deriv = 2. * Aa_;
            break;
        case PhaseIndex::phaseB:
            deriv = 2. * Ab_;
            break;
        default:
            return;
    }

    d2fdc2[0] = deriv;
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature

bool QuadraticFreeEnergyFunctionsBinaryThreePhase::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
    return false;
}

//=======================================================================

void QuadraticFreeEnergyFunctionsBinaryThreePhase::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc,
    double& fl, double& fa, double& fb)
{
    // std::cout<<"QuadraticFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;

    double c[3] = { conc, conc, conc };

    QuadraticConcSolverBinaryThreePhase solver;
    solver.setup(
        conc, hphi[0], hphi[1], hphi[2], Al_, ceql_, Aa_, ceqa_, Ab_, ceqb_);
    int ret = solver.ComputeConcentration(c);
    if (ret < 0)
    {
#if 0
        std::cerr << "ERROR in "
                     "QuadraticFreeEnergyFunctionsBinaryThreePhase::"
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

int QuadraticFreeEnergyFunctionsBinaryThreePhase::computePhaseConcentrations(
    const double temperature, const double* const conc, const double* const phi,
    double* x)
{
    const double hphi0 = interp_func(conc_interp_func_type_, phi[0]);
    const double hphi1 = interp_func(conc_interp_func_type_, phi[1]);
    const double hphi2 = interp_func(conc_interp_func_type_, phi[2]);

    // solve system of equations to find (cl,cs) given conc[0] and hphi
    // x: initial guess and solution
    QuadraticConcSolverBinaryThreePhase solver;
    solver.setup(
        conc[0], hphi0, hphi1, hphi2, Al_, ceql_, Aa_, ceqa_, Ab_, ceqb_);
    int ret = solver.ComputeConcentration(x);
#if 0
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "QuadraticFreeEnergyFunctionsBinary::"
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

//-----------------------------------------------------------------------
void QuadraticFreeEnergyFunctionsBinaryThreePhase::energyVsPhiAndC(
    const double temperature, const double* const ceq, const bool found_ceq,
    const double phi_well_scale, const int npts_phi, const int npts_c)
{
    // Not implemented because it is ill-defined for a three-phase system.
}

// Print out free energy as a function of phase
// for given composition and temperature
// File format: ASCII VTK, readble with Visit
void QuadraticFreeEnergyFunctionsBinaryThreePhase::printEnergyVsPhiHeader(
    const double temperature, const int nphi, const int nc, const double cmin,
    const double cmax, const double slopec, std::ostream& os) const
{
    // Not implemented because it is ill-defined for a three-phase system
}

//=======================================================================
void QuadraticFreeEnergyFunctionsBinaryThreePhase::printEnergyVsPhi(
    const double* const conc, const double temperature,
    const double phi_well_scale, const int npts, const double slopec,
    std::ostream& os)
{
    // Not implemented because it is ill-defined for a three-phase system
}

//=======================================================================
// compute free energy in [J/mol]

double QuadraticFreeEnergyFunctionsBinaryThreePhase::fchem(
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

void QuadraticFreeEnergyFunctionsBinaryThreePhase::printEnergyVsComposition(
    const double temperature, std::ostream& os, const double cmin,
    const double cmax, const int npts)
{
    const double dc = (cmax - cmin) / (double)(npts - 1);

    os << "#phi0=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        const double phi[3] = { 1., 0., 0. };

        double e = fchem(phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi1=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        const double phi[3] = { 0., 1., 0. };
        double e            = fchem(phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi2=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        const double phi[3] = { 0., 0., 1. };
        double e            = fchem(phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
}

//=======================================================================

void QuadraticFreeEnergyFunctionsBinaryThreePhase::preRunDiagnostics(
    const double T0, const double T1)
{
}
}
