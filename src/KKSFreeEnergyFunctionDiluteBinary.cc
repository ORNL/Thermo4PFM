#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "PhysicalConstants.h"
#include "functions.h"
#include "well_functions.h"
#include "xlogx.h"

#include <iostream>
#include <math.h>

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

KKSFreeEnergyFunctionDiluteBinary::KKSFreeEnergyFunctionDiluteBinary(
    pt::ptree& conc_db, const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type)
    : energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type),
      tol_(1.e-8),
      maxiters_(20),
      alpha_(1.)
{
    std::string fenergy_diag_filename("energy.vtk");
    fenergy_diag_filename_ = new char[fenergy_diag_filename.length() + 1];
    strcpy(fenergy_diag_filename_, fenergy_diag_filename.c_str());

    readParameters(conc_db);

    fA_ = log(1. / ke_);

    boost::optional<pt::ptree&> newton_db;

    if (newton_db) readNewtonparameters(newton_db.get());
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::readNewtonparameters(
    pt::ptree& newton_db)
{
    tol_      = newton_db.get<double>("tol", tol_);
    alpha_    = newton_db.get<double>("alpha", alpha_);
    maxiters_ = newton_db.get<int>("max_its", maxiters_);
    assert(maxiters_ > 1);
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::readParameters(pt::ptree& conc_db)
{
    me_ = conc_db.get<double>("liquidus_slope");
    Tm_ = conc_db.get<double>("meltingT");
    ke_ = conc_db.get<double>("keq");

    assert(me_ < 0.);
    assert(Tm_ > 0.);
    assert(ke_ <= 1.);
}

//-----------------------------------------------------------------------

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

double KKSFreeEnergyFunctionDiluteBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    double fe = xlogx(conc[0]) + xlogx(1. - conc[0]);
    double fB = computeFB(temperature);

    switch (pi)
    {
        case PhaseIndex::phaseL:
            break;
        case PhaseIndex::phaseA:
            fe += conc[0] * fA_ + (1. - conc[0]) * fB;
            break;
        default:
#ifndef HAVE_OPENMP_OFFLOAD
            std::cout << "KKSFreeEnergyFunctionDiluteBinary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
            abort();
#endif
            return 0.;
    }

    fe *= GASCONSTANT_R_JPKPMOL * temperature;

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

void KKSFreeEnergyFunctionDiluteBinary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    double mu = xlogx_deriv(conc[0]) - xlogx_deriv(1.0 - conc[0]);

    switch (pi)
    {
        case PhaseIndex::phaseL:
            break;
        case PhaseIndex::phaseA:
            mu += (fA_ - computeFB(temperature));
            break;
        default:
#ifndef HAVE_OPENMP_OFFLOAD
            std::cout << "KKSFreeEnergyFunctionDiluteBinary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
            abort();
#endif
            return;
    }

    deriv[0] = GASCONSTANT_R_JPKPMOL * temperature * mu;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    double* d2fdc2)
{
    (void)pi;

#ifndef HAVE_OPENMP_OFFLOAD
    assert(conc[0] >= 0.);
    assert(conc[0] <= 1.);
#endif

    const double rt = GASCONSTANT_R_JPKPMOL * temp;

    d2fdc2[0] = rt * (xlogx_deriv2(conc[0]) + xlogx_deriv2(1.0 - conc[0]));
}

//=======================================================================

double KKSFreeEnergyFunctionDiluteBinary::computeFB(
    const double temperature) const
{
#ifndef HAVE_OPENMP_OFFLOAD
    assert(ke_ > 0.);
    assert(temperature < Tm_);
    assert(me_ < 0.);
#endif

    const double cLe = (temperature - Tm_) / me_;
    const double cSe = cLe * ke_;

#ifndef HAVE_OPENMP_OFFLOAD
    assert(cLe < 1.);
    assert(cSe < 1.);
#endif

    return log(1. - cLe) - log(1. - cSe);
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool KKSFreeEnergyFunctionDiluteBinary::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
    (void)maxits;

#ifndef HAVE_OPENMP_OFFLOAD
    if (verbose)
        std::cout << "KKSFreeEnergyFunctionDiluteBinary::computeCeqT()"
                  << std::endl;
    assert(temperature > 0.);
#endif

    ceq[0] = (temperature - Tm_) / me_;
    ceq[1] = ceq[0] * ke_;

    return true;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc,
    double& fl, double& fa)
{
    // std::cout<<"KKSFreeEnergyFunctionDiluteBinary::computePhasesFreeEnergies()"<<endl;

    double c[2] = { conc, conc };

    double fB = computeFB(temperature);

    KKSdiluteBinaryConcSolver solver;
    solver.setup(conc, *hphi, 1. - *hphi, fA_, fB);
    int ret = solver.ComputeConcentration(c, tol_, maxiters_, alpha_);

#ifndef HAVE_OPENMP_OFFLOAD
    if (ret < 0)
    {
        std::cerr << "ERROR in "
                     "KKSFreeEnergyFunctionDiluteBinary::"
                     "computePhasesFreeEnergies() "
                     "---"
                  << "conc=" << conc << ", hphi=" << &hphi << std::endl;
        abort();
    }
    assert(c[0] >= 0.);
    assert(c[1] >= 0.);
#endif

    fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);
    fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------

int KKSFreeEnergyFunctionDiluteBinary::computePhaseConcentrations(
    const double temperature, const double* const conc,
    const double* const hphi, double* x)

{
#ifndef HAVE_OPENMP_OFFLOAD
    assert(x[0] >= 0.);
    assert(x[1] >= 0.);
    assert(x[0] <= 1.);
    assert(x[1] <= 1.);
    assert(maxiters_ > 1);
#endif

    const double fB = computeFB(temperature);

    // conc could be outside of [0.,1.] in a trial step
    double c0 = conc[0] >= 0. ? conc[0] : 0.;
    c0        = c0 <= 1. ? c0 : 1.;
    KKSdiluteBinaryConcSolver solver;
    solver.setup(c0, hphi[0], 1. - hphi[0], fA_, fB);
    int ret = solver.ComputeConcentration(x, tol_, maxiters_);

    return ret;
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

//-----------------------------------------------------------------------

void KKSFreeEnergyFunctionDiluteBinary::energyVsPhiAndC(
    const double temperature, const double* const ceq, const bool found_ceq,
    const double phi_well_scale, const int npts_phi, const int npts_c)
{
    std::cout << "KKSFreeEnergyFunctionDiluteBinary::energyVsPhiAndC()..."
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
    std::cout << "KKSFreeEnergyFunctionDiluteBinary: Use slope: " << slopec
              << "..." << std::endl;

    {
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
}

// Print out free energy as a function of phase
// for given composition and temperature
// File format: ASCII VTK, readble with Visit
void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsPhiHeader(
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

void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsPhi(
    const double* const conc, const double temperature,
    const double phi_well_scale, const int npts, const double slopec,
    std::ostream& os)
{
    // std::cout << "KKSFreeEnergyFunctionDiluteBinary::printEnergyVsPhi()..."
    // << std::endl;
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
double KKSFreeEnergyFunctionDiluteBinary::fchem(
    const double* const phi, const double* const conc, const double temperature)
{
    const double hcphi = interp_func(conc_interp_func_type_, phi[0]);

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;
    if ((phi[0] > tol) & (phi[0] < (1. - tol)))
    {
        computePhasesFreeEnergies(temperature, &hcphi, conc[0], fl, fa);
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
    double e           = (1.0 - hfphi) * fl + hfphi * fa;

    return e;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsComposition(
    const double temperature, std::ostream& os, const int npts)
{
    const double dc = 1.0 / (double)(npts - 1);

    os << "#phi=0" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        const double phi = 0.0;

        double e = fchem(&phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        const double phi = 1.0;

        double e = fchem(&phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
}
}
