#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADConcSolverBinary.h"
#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADFunctions.h"
#include "PhysicalConstants.h"
#include "functions.h"
#include "well_functions.h"

#include <boost/property_tree/json_parser.hpp>

#include <iomanip>
#include <string>

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

CALPHADFreeEnergyFunctionsBinary::CALPHADFreeEnergyFunctionsBinary(
    pt::ptree& calphad_db, boost::optional<pt::ptree&> newton_db,
    const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type)
    : energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type),
      newton_tol_(1.e-8),
      newton_alpha_(1.),
      newton_maxits_(20),
      newton_verbose_(false)
{
    std::string fenergy_diag_filename("energy.vtk");
    fenergy_diag_filename_ = new char[fenergy_diag_filename.length() + 1];
    strcpy(fenergy_diag_filename_, fenergy_diag_filename.c_str());

    readParameters(calphad_db);

    if (newton_db) readNewtonparameters(newton_db.get());
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::readNewtonparameters(
    pt::ptree& newton_db)
{
    newton_tol_     = newton_db.get<double>("tol", newton_tol_);
    newton_alpha_   = newton_db.get<double>("alpha", newton_alpha_);
    newton_maxits_  = newton_db.get<int>("max_its", newton_maxits_);
    newton_verbose_ = newton_db.get<bool>("verbose", newton_verbose_);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::readParameters(pt::ptree& calphad_db)
{
    pt::ptree& species0_db = calphad_db.get_child("SpeciesA");
    std::string dbnameL("PhaseL");
    g_species_phaseL_[0].initialize("L0", species0_db.get_child(dbnameL));
    std::string dbnameA("PhaseA");
    g_species_phaseA_[0].initialize("A0", species0_db.get_child(dbnameA));

    pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
    g_species_phaseL_[1].initialize("L1", speciesB_db.get_child(dbnameL));
    g_species_phaseA_[1].initialize("A1", speciesB_db.get_child(dbnameA));

    // read Lmix coefficients
    std::string dbnamemixL("LmixPhaseL");
    pt::ptree Lmix0_db = calphad_db.get_child(dbnamemixL);
    readLmixBinary(Lmix0_db, LmixPhaseL_);

    std::string dbnamemixA("LmixPhaseA");
    pt::ptree Lmix1_db = calphad_db.get_child(dbnamemixA);
    readLmixBinary(Lmix1_db, LmixPhaseA_);

    // print database just read
    // std::clog << "CALPHAD database..." << std::endl;
    // pt::write_json(std::clog, calphad_db);
}

//-----------------------------------------------------------------------
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

double CALPHADFreeEnergyFunctionsBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    const double l0 = lmixPhase(0, pi, temperature);
    const double l1 = lmixPhase(1, pi, temperature);
    const double l2 = lmixPhase(2, pi, temperature);
    const double l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &g_species_phaseL_[0];
            break;
        case PhaseIndex::phaseA:
            g_species = &g_species_phaseA_[0];
            break;
        default:
            //            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
            //                         "computeFreeEnergy(), undefined phase"
            //                      << "!!!" << std::endl;
            abort();
            return 0.;
    }

    double fe = conc[0] * g_species[0].fenergy(temperature)
                + (1. - conc[0]) * g_species[1].fenergy(temperature)
                + CALPHADcomputeFMixBinary(l0, l1, l2, l3, conc[0])
                + CALPHADcomputeFIdealMixBinary(
                      gas_constant_R_JpKpmol * temperature, conc[0]);

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

void CALPHADFreeEnergyFunctionsBinary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    const double l0 = lmixPhase(0, pi, temperature);
    const double l1 = lmixPhase(1, pi, temperature);
    const double l2 = lmixPhase(2, pi, temperature);
    const double l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &g_species_phaseL_[0];
            break;
        case PhaseIndex::phaseA:
            g_species = &g_species_phaseA_[0];
            break;
        default:
            //            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
            //                         "computeFreeEnergy(), undefined phase!!!"
            //                      << std::endl;
            abort();
            return;
    }

    double mu = (g_species[0].fenergy(temperature)
                    - g_species[1].fenergy(temperature))
                + CALPHADcomputeFMix_derivBinary(l0, l1, l2, l3, conc[0])
                + CALPHADcomputeFIdealMix_derivBinary(
                      gas_constant_R_JpKpmol * temperature, conc[0]);

    deriv[0] = mu;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    double* d2fdc2)
{
    // assert(conc[0] >= 0.);
    // assert(conc[0] <= 1.);

    const double l0 = lmixPhase(0, pi, temp);
    const double l1 = lmixPhase(1, pi, temp);
    const double l2 = lmixPhase(2, pi, temp);
    const double l3 = lmixPhase(3, pi, temp);
    const double rt = gas_constant_R_JpKpmol * temp;

    d2fdc2[0] = (CALPHADcomputeFMix_deriv2Binary(l0, l1, l2, l3, conc[0])
                 + CALPHADcomputeFIdealMix_deriv2Binary(rt, conc[0]));
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computeTdependentParameters(
    const double temperature, double* Lmix_L, double* Lmix_A, double* fA,
    double* fB)
{
    fA[0]     = g_species_phaseL_[0].fenergy(temperature);
    fB[0]     = g_species_phaseL_[1].fenergy(temperature);
    Lmix_L[0] = lmixPhase(0, PhaseIndex::phaseL, temperature);
    Lmix_L[1] = lmixPhase(1, PhaseIndex::phaseL, temperature);
    Lmix_L[2] = lmixPhase(2, PhaseIndex::phaseL, temperature);
    Lmix_L[3] = lmixPhase(3, PhaseIndex::phaseL, temperature);

    fA[1]     = g_species_phaseA_[0].fenergy(temperature);
    fB[1]     = g_species_phaseA_[1].fenergy(temperature);
    Lmix_A[0] = lmixPhase(0, PhaseIndex::phaseA, temperature);
    Lmix_A[1] = lmixPhase(1, PhaseIndex::phaseA, temperature);
    Lmix_A[2] = lmixPhase(2, PhaseIndex::phaseA, temperature);
    Lmix_A[3] = lmixPhase(3, PhaseIndex::phaseA, temperature);
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyFunctionsBinary::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
#ifndef HAVE_OPENMP_OFFLOAD
    if (verbose)
        std::cout << "CALPHADFreeEnergyFunctionsBinary::computeCeqT()"
                  << std::endl;
#endif
    // assert(temperature > 0.);

    // evaluate temperature dependent parameters
    double fA[3];
    double fB[3];

    double Lmix_L[4];
    double Lmix_A[4];

    computeTdependentParameters(temperature, Lmix_L, Lmix_A, fA, fB);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    CALPHADEqConcSolverBinary eq_solver;

    eq_solver.setup(RTinv, Lmix_L, Lmix_A, fA, fB);
    int ret = eq_solver.ComputeConcentration(ceq, newton_tol_, maxits);

#ifndef HAVE_OPENMP_OFFLOAD
    if (ret >= 0)
    {
        if (verbose)
        {
            std::cout << "CALPHAD, c_eq phase0=" << ceq[0] << std::endl;
            std::cout << "CALPHAD, c_eq phase1=" << ceq[1] << std::endl;
        }
    }
    else
    {
        std::cout << "CALPHADFreeEnergyFunctionsBinary, WARNING: ceq "
                     "computation did not converge"
                  << std::endl;
    }
#endif

    return (ret >= 0);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies(
    const double temperature, const double hphi, const double conc, double& fl,
    double& fa)
{
    // std::cout<<"CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;

    double c[2] = { conc, conc };

    // evaluate temperature dependent parameters
    double fA[2];
    double fB[2];

    double Lmix_L[4];
    double Lmix_A[4];
    computeTdependentParameters(temperature, Lmix_L, Lmix_A, fA, fB);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    CALPHADConcSolverBinary solver;
    solver.setup(conc, hphi, RTinv, Lmix_L, Lmix_A, fA, fB);
    int ret = solver.ComputeConcentration(
        c, newton_tol_, newton_maxits_, newton_alpha_);
    if (ret < 0)
    {
#if 0
        std::cerr << "ERROR in "
                     "CALPHADFreeEnergyFunctionsBinary::"
                     "computePhasesFreeEnergies()"
                     " ---"
                  << "conc=" << conc << ", hphi=" << hphi << std::endl;
#endif
        abort();
    }

    // assert(c[0] >= 0.);
    fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);

    // assert(c[1] >= 0.);
    fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------

int CALPHADFreeEnergyFunctionsBinary::computePhaseConcentrations(
    const double temperature, const double* const conc, const double phi,
    double* x)
{
    // assert(x[0] >= 0.);
    // assert(x[1] >= 0.);
    // assert(x[0] <= 1.);
    // assert(x[1] <= 1.);

    const double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    double fA[2];
    double fB[2];
    double Lmix_L[4];
    double Lmix_A[4];

    computeTdependentParameters(temperature, Lmix_L, Lmix_A, fA, fB);

    const double hphi = interp_func(conc_interp_func_type_, phi);

    // conc could be outside of [0.,1.] in a trial step
    double c0 = conc[0] >= 0. ? conc[0] : 0.;
    c0        = c0 <= 1. ? c0 : 1.;
    // solve system of equations to find (cl,cs) given c0 and hphi
    // x: initial guess and solution
    CALPHADConcSolverBinary solver;
    solver.setup(c0, hphi, RTinv, Lmix_L, Lmix_A, fA, fB);
    int ret = solver.ComputeConcentration(
        x, newton_tol_, newton_maxits_, newton_alpha_);
#if 0
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "CALPHADFreeEnergyFunctionsBinary::"
                     "computePhaseConcentrations() "
                     "failed for conc="
                  << conc[0] << ", hphi=" << hphi << std::endl;
        abort();
    }
#endif

    return ret;
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

//-----------------------------------------------------------------------

void CALPHADFreeEnergyFunctionsBinary::energyVsPhiAndC(const double temperature,
    const double* const ceq, const bool found_ceq, const double phi_well_scale,
    const int npts_phi, const int npts_c)
{
    std::cout << "CALPHADFreeEnergyFunctionsBinary::energyVsPhiAndC()..."
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
    std::cout << "CALPHADFreeEnergyFunctionsBinary: Use slope: " << slopec
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
void CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhiHeader(
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

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhi(
    const double* const conc, const double temperature,
    const double phi_well_scale, const int npts, const double slopec,
    std::ostream& os)
{
    // std::cout << "CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhi()..." <<
    // std::endl;
    const double dphi = 1.0 / (double)(npts - 1);

    // os << "# phi     f(phi)     for c=" << conc
    //           << " and T=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double phi = i * dphi;

        double e       = fchem(phi, conc, temperature);
        const double w = phi_well_scale * well_func(phi);

        os << e + w + slopec * conc[0] << std::endl;
    }
    // os << std::endl;
}

//=======================================================================
// compute free energy in [J/mol]
double CALPHADFreeEnergyFunctionsBinary::fchem(
    const double phi, const double* const conc, const double temperature)
{
    const double hcphi = interp_func(conc_interp_func_type_, phi);

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;
    double fb        = 0.;
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

    const double hfphi = interp_func(energy_interp_func_type_, phi);

    double e = (1.0 - hfphi) * fl + hfphi * fa;

    return e;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsComposition(
    const double temperature, std::ostream& os, const int npts)
{
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

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::preRunDiagnostics(
    const double T0, const double T1)
{
    std::ofstream os1("FlC0vsT.dat", std::ios::out);
    os1 << "#Species 0, Phase L" << std::endl;
    g_species_phaseL_[0].plotFofT(os1, T0, T1);

    std::ofstream os2("FlC1vsT.dat", std::ios::out);
    os2 << "#Species 1, Phase L" << std::endl;
    g_species_phaseL_[1].plotFofT(os2, T0, T1);

    std::ofstream os3("FsC0vsT.dat", std::ios::out);
    os3 << "#Species 0, Phase A" << std::endl;
    g_species_phaseA_[0].plotFofT(os3, T0, T1);

    std::ofstream os4("FsC1vsT.dat", std::ios::out);
    os4 << "#Species 1, Phase A" << std::endl;
    g_species_phaseA_[1].plotFofT(os4, T0, T1);
}
}
