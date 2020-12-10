#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "PhysicalConstants.h"
#include "functions.h"
#include "well_functions.h"
#include "xlogx.h"

#include <cmath>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

KKSFreeEnergyFunctionDiluteBinary::KKSFreeEnergyFunctionDiluteBinary(
    pt::ptree& conc_db, const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type)
    : energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type)
{
    fenergy_diag_filename_ = "energy.vtk";

    ceq_l_ = -1;
    ceq_a_ = -1;

    readParameters(conc_db);

    fA_ = log(1. / ke_);

    boost::optional<pt::ptree&> newton_db;

    setupSolver(newton_db);
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::setupSolver(
    boost::optional<pt::ptree&> newton_db)
{
    std::cout << "KKSFreeEnergyFunctionDiluteBinary::setupSolver()..."
              << std::endl;
    solver_ = new KKSdiluteBinaryConcentrationSolver();

    if (newton_db) readNewtonparameters(newton_db.get());
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::readNewtonparameters(
    pt::ptree& newton_db)
{
    double tol         = newton_db.get<double>("tol", 1.e-8);
    double alpha       = newton_db.get<double>("alpha", 1.);
    int maxits         = newton_db.get<int>("max_its", 20);
    const bool verbose = newton_db.get<bool>("verbose", false);

    solver_->SetTolerance(tol);
    solver_->SetMaxIterations(maxits);
    solver_->SetDamping(alpha);
    solver_->SetVerbose(verbose);
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

double KKSFreeEnergyFunctionDiluteBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    double fe = xlogx(conc[0]) + xlogx(1. - conc[0]);
    setupFB(temperature);

    switch (pi)
    {
        case PhaseIndex::phaseL:
            break;
        case PhaseIndex::phaseA:
            fe += conc[0] * fA_ + (1. - conc[0]) * fB_;
            break;
        default:
            std::cout << "KKSFreeEnergyFunctionDiluteBinary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
            abort();
            return 0.;
    }

    fe *= gas_constant_R_JpKpmol * temperature;

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
            setupFB(temperature);
            mu += (fA_ - fB_);
            break;
        default:
            std::cout << "KKSFreeEnergyFunctionDiluteBinary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
            abort();
            return;
    }

    deriv[0] = gas_constant_R_JpKpmol * temperature * mu;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    std::vector<double>& d2fdc2)
{
    assert(conc[0] >= 0.);
    assert(conc[0] <= 1.);

    const double rt = gas_constant_R_JpKpmol * temp;

    d2fdc2[0] = rt * (xlogx_deriv2(conc[0]) + xlogx_deriv2(1.0 - conc[0]));
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::setupFB(const double temperature)
{
    assert(ke_ > 0.);
    assert(temperature < Tm_);
    assert(me_ < 0.);

    const double cLe = (temperature - Tm_) / me_;
    const double cSe = cLe * ke_;

    assert(cLe < 1.);
    assert(cSe < 1.);

    fB_ = std::log(1. - cLe) - std::log(1. - cSe);
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool KKSFreeEnergyFunctionDiluteBinary::computeCeqT(const double temperature,
    const PhaseIndex pi0, const PhaseIndex pi1, double* ceq, const int maxits,
    const bool verbose)
{
    if (verbose)
        std::cout << "KKSFreeEnergyFunctionDiluteBinary::computeCeqT()"
                  << std::endl;
    assert(temperature > 0.);

    ceq_l_ = (temperature - Tm_) / me_;
    ceq_a_ = ceq_l_ * ke_;

    ceq[0] = ceq_l_;
    ceq[1] = ceq_a_;

    return true;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::computePhasesFreeEnergies(
    const double temperature, const double hphi, const double conc, double& fl,
    double& fa)
{
    // std::cout<<"KKSFreeEnergyFunctionDiluteBinary::computePhasesFreeEnergies()"<<endl;

    double c[2] = { conc, conc };
    if (ceq_l_ >= 0.) c[0] = ceq_l_;
    if (ceq_a_ >= 0.) c[1] = ceq_a_;

    setupFB(temperature);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);
    int ret = solver_->ComputeConcentration(c, conc, hphi, RTinv, fA_, fB_);

    if (ret < 0)
    {
        std::cerr << "ERROR in "
                     "KKSFreeEnergyFunctionDiluteBinary::"
                     "computePhasesFreeEnergies() "
                     "---"
                  << "conc=" << conc << ", hphi=" << hphi << std::endl;
        abort();
    }

    assert(c[0] >= 0.);
    fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);

    assert(c[1] >= 0.);
    fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------

int KKSFreeEnergyFunctionDiluteBinary::computePhaseConcentrations(
    const double temperature, const double* const conc, const double phi,
    double* x)

{
    assert(x[0] >= 0.);
    assert(x[1] >= 0.);
    assert(x[0] <= 1.);
    assert(x[1] <= 1.);

    const double conc0 = conc[0];

    const char interp_func_type = concInterpChar(conc_interp_func_type_);
    const double hphi           = interp_func(phi, interp_func_type);

    setupFB(temperature);

    // conc could be outside of [0.,1.] in a trial step
    double c0 = conc[0] >= 0. ? conc[0] : 0.;
    c0        = c0 <= 1. ? c0 : 1.;
    int ret   = solver_->ComputeConcentration(x, c0, hphi,
        -1., // unused parameter
        fA_, fB_);
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "KKSFreeEnergyFunctionDiluteBinary::"
                     "computePhaseConcentrations() "
                     "failed for conc="
                  << conc0 << ", hphi=" << hphi << std::endl;
        abort();
    }

    return ret;
}

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

        std::ofstream tfile(fenergy_diag_filename_.data(), std::ios::out);

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

        double e       = fchem(phi, conc, temperature);
        const double w = phi_well_scale * well_func(phi);

        os << e + w + slopec * conc[0] << std::endl;
    }
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsEta(
    const double* const conc, const double temperature, const int npts,
    const double slopec, std::ostream& os)
{
    (void)conc;
    (void)temperature;
    (void)npts;
    (void)slopec;
    (void)os;
}

//=======================================================================
// compute free energy in [J/mol]
double KKSFreeEnergyFunctionDiluteBinary::fchem(
    const double phi, const double* const conc, const double temperature)
{

    const char interp_func_type = concInterpChar(conc_interp_func_type_);
    const double hcphi          = interp_func(phi, interp_func_type);

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;
    if ((phi > tol) & (phi < (1. - tol)))
    {
        computePhasesFreeEnergies(temperature, hcphi, conc[0], fl, fa);
    }
    else
    {
        if (phi <= tol)
        {
            fl = computeFreeEnergy(temperature, conc, PhaseIndex::phaseL);
        }
        else
        {
            fa = computeFreeEnergy(temperature, conc, PhaseIndex::phaseA);
        }
    }

    const char interpf = energyInterpChar(energy_interp_func_type_);
    const double hfphi = interp_func(phi, interpf);
    double e           = (1.0 - hfphi) * fl + hfphi * fa;

    return e;
}

//=======================================================================

void KKSFreeEnergyFunctionDiluteBinary::printEnergyVsComposition(
    const double temperature, const int npts)
{
    std::ofstream os("FvsC.dat", std::ios::out);

    const double dc = 1.0 / (double)(npts - 1);

    os << "#phi=0" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        double e = fchem(0., &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        double e = fchem(1., &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
}
}
