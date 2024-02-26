#include "QuadraticFreeEnergyFunctionsBinary.h"
#include "functions.h"
#include "well_functions.h"

#include <iomanip>
#include <string>

namespace Thermo4PFM
{

QuadraticFreeEnergyFunctionsBinary::QuadraticFreeEnergyFunctionsBinary(
    const double Tref, const double Al, const double ceql,
    const double m_liquid, const double Aa, const double ceqa,
    const double m_solid, const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type)
    : Tref_(Tref),
      Al_(Al),
      ceql_(ceql),
      m_liquid_(m_liquid),
      Aa_(Aa),
      ceqa_(ceqa),
      m_solid_(m_solid),
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
double QuadraticFreeEnergyFunctionsBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    double A;
    double ceq;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            A   = Al_;
            ceq = ceql_ + (temperature - Tref_) * m_liquid_;
            break;
        case PhaseIndex::phaseA:
            A   = Aa_;
            ceq = ceqa_ + (temperature - Tref_) * m_solid_;
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

void QuadraticFreeEnergyFunctionsBinary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    double A;
    double ceq;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            A   = Al_;
            ceq = ceql_ + (temperature - Tref_) * m_liquid_;
            break;
        case PhaseIndex::phaseA:
            A   = Aa_;
            ceq = ceqa_ + (temperature - Tref_) * m_solid_;
            break;
        default:
            return;
    }

    double mu = 2. * A * (conc[0] - ceq);

    deriv[0] = mu;
}

//=======================================================================

void QuadraticFreeEnergyFunctionsBinary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    double* d2fdc2)
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
        default:
            return;
    }

    d2fdc2[0] = deriv;
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature

bool QuadraticFreeEnergyFunctionsBinary::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
    ceq[0] = ceql_ + (temperature - Tref_) * m_liquid_;
    ceq[1] = ceqa_ + (temperature - Tref_) * m_solid_;

    return true;
}

//=======================================================================

double QuadraticFreeEnergyFunctionsBinary::computeLiquidConcentration(
    const double temp, const double hphi, const double c) const
{
    double Ceq_liquid  = ceql_ + (temp - Tref_) * m_liquid_;
    double Ceq_solid_A = ceqa_ + (temp - Tref_) * m_solid_;

    return (c - hphi * (Ceq_solid_A - (Al_ / Aa_) * Ceq_liquid))
           / ((1.0 - hphi) + hphi * (Al_ / Aa_));
}

//=======================================================================

double QuadraticFreeEnergyFunctionsBinary::computeSolidAConcentration(
    const double temp, const double hphi, const double c) const
{
    double Ceq_liquid  = ceql_ + (temp - Tref_) * m_liquid_;
    double Ceq_solid_A = ceqa_ + (temp - Tref_) * m_solid_;

    return (c - (1.0 - hphi) * (Ceq_liquid - (Aa_ / Al_) * Ceq_solid_A))
           / ((1.0 - hphi) * (Aa_ / Al_) + hphi);
}

//=======================================================================

void QuadraticFreeEnergyFunctionsBinary::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc,
    double& fl, double& fa)
{
    double c[1] = { conc };
    double x[2];
    computePhaseConcentrations(temperature, c, hphi, x);

    fl = computeFreeEnergy(temperature, &x[0], PhaseIndex::phaseL, false);
    fa = computeFreeEnergy(temperature, &x[1], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------

int QuadraticFreeEnergyFunctionsBinary::computePhaseConcentrations(
    const double temperature, const double* const conc, const double* const phi,
    double* x)
{
    // solve system of equations to find (cl,cs) given conc[0] and hphi
    // x: initial guess and solution
    x[0] = computeLiquidConcentration(temperature, phi[0], conc[0]);
    x[1] = computeSolidAConcentration(temperature, phi[0], conc[0]);

    return 0;
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

//-----------------------------------------------------------------------
void QuadraticFreeEnergyFunctionsBinary::energyVsPhiAndC(
    const double temperature, const double* const ceq, const bool found_ceq,
    const double phi_well_scale, const int npts_phi, const int npts_c)
{
    std::cout << "QuadraticFreeEnergyFunctionsBinary::energyVsPhiAndC()..."
              << std::endl;

    double slopec = 0.;
    double fc0    = 0.;
    double fc1    = 0.;
    if (found_ceq)
    {
        // compute slope of f between equilibrium concentrations
        // to add slopec*conc to energy later on

        fc0    = computeFreeEnergy(temperature, &ceq[0], PhaseIndex::phaseL);
        fc1    = computeFreeEnergy(temperature, &ceq[1], PhaseIndex::phaseA);
        slopec = -(fc1 - fc0) / (ceq[1] - ceq[0]);
    }
    std::cout << std::setprecision(8) << "fc0: " << fc0 << "..."
              << ", fc1: " << fc1 << "..." << std::endl;
    std::cout << "QuadraticFreeEnergyFunctionsBinary: Use slope: " << slopec
              << "..." << std::endl;

    // reset cmin, cmax, deltac
    double cmin   = std::min(ceq[0], ceq[1]);
    double cmax   = std::max(ceq[0], ceq[1]);
    double dc     = cmax - cmin;
    cmin          = std::max(0.25 * cmin, cmin - 0.25 * dc);
    cmax          = std::min(1. - 0.25 * (1. - cmax), cmax + 0.25 * dc);
    cmax          = std::max(cmax, cmin + dc);
    double deltac = (cmax - cmin) / (npts_c - 1);

    std::ofstream tfile(fenergy_diag_filename_, std::ios::out);

    printEnergyVsPhiHeader(
        temperature, npts_phi, npts_c, cmin, cmax, slopec, tfile);

    for (int i = 0; i < npts_c; i++)
    {
        double conc = cmin + deltac * i;
        printEnergyVsPhi(
            &conc, temperature, phi_well_scale, npts_phi, slopec, tfile);
    }
}

// Print out free energy as a function of phase
// for given composition and temperature
// File format: ASCII VTK, readble with Visit
void QuadraticFreeEnergyFunctionsBinary::printEnergyVsPhiHeader(
    const double temperature, const int nphi, const int nc, const double cmin,
    const double cmax, const double slopec, std::ostream& os) const
{
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << "Free energy + " << slopec << "*c [J/mol] at T=" << temperature
       << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET STRUCTURED_POINTS" << std::endl;

    os << "DIMENSIONS   " << nphi << " " << nc << " 1" << std::endl;
    double asp_ratio_c = (nc > 1) ? (cmax - cmin) / (nc - 1) : 1.;
    os << "ASPECT_RATIO " << 1. / (nphi - 1) << " " << asp_ratio_c << " 1."
       << std::endl;
    os << "ORIGIN        0. " << cmin << " 0." << std::endl;
    os << "POINT_DATA   " << nphi * nc << std::endl;
    os << "SCALARS energy float 1" << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
}

//=======================================================================
void QuadraticFreeEnergyFunctionsBinary::printEnergyVsPhi(
    const double* const conc, const double temperature,
    const double phi_well_scale, const int npts, const double slopec,
    std::ostream& os)
{
    const double dphi = 1.0 / (double)(npts - 1);

    for (int i = 0; i < npts; i++)
    {
        const double phi = i * dphi;

        double e       = fchem(&phi, conc, temperature);
        const double w = phi_well_scale * well_func(phi);

        os << e + w + slopec * conc[0] << std::endl;
    }
}

//=======================================================================
// compute free energy in [J/mol]

double QuadraticFreeEnergyFunctionsBinary::fchem(
    const double* const phi, const double* const conc, const double temperature)
{
    double hcphi[1];
    hcphi[0] = interp_func(conc_interp_func_type_, phi[0]);

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;

    if ((1.0 - phi[0] > tol) && (1.0 - phi[1] > tol) && (1.0 - phi[2] > tol))
    {
        computePhasesFreeEnergies(temperature, hcphi, conc[0], fl, fa);
    }
    else
    {
        if (phi[0] <= tol)
        {
            fl = computeFreeEnergy(temperature, conc, PhaseIndex::phaseL);
        }
        else
        {
            fa = computeFreeEnergy(temperature, conc, PhaseIndex::phaseA);
        }
    }

    const double hfphi = interp_func(energy_interp_func_type_, phi[0]);

    return (1.0 - hfphi) * fl + hfphi * fa;
}

//=======================================================================

void QuadraticFreeEnergyFunctionsBinary::printEnergyVsComposition(
    const double temperature, std::ostream& os, const double cmin,
    const double cmax, const int npts)
{
    const double dc = (cmax - cmin) / (double)(npts - 1);

    os << "#phi0=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        const double phi = 0.0;

        double e = fchem(&phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi1=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        const double phi = 1.0;

        double e = fchem(&phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;
}

//=======================================================================

void QuadraticFreeEnergyFunctionsBinary::preRunDiagnostics(
    const double T0, const double T1)
{
}
}
