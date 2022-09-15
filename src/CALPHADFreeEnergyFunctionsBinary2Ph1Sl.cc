#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "CALPHADConcSolverBinary2Ph1Sl.h"
#include "CALPHADFunctions.h"
#include "PhysicalConstants.h"
#include "functions.h"
#include "well_functions.h"

#include <boost/property_tree/json_parser.hpp>

#include <iomanip>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

CALPHADFreeEnergyFunctionsBinary2Ph1Sl::CALPHADFreeEnergyFunctionsBinary2Ph1Sl(
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

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::readNewtonparameters(
    pt::ptree& newton_db)
{
    newton_tol_     = newton_db.get<double>("tol", newton_tol_);
    newton_alpha_   = newton_db.get<double>("alpha", newton_alpha_);
    newton_maxits_  = newton_db.get<int>("max_its", newton_maxits_);
    newton_verbose_ = newton_db.get<bool>("verbose", newton_verbose_);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::readParameters(
    pt::ptree& calphad_db)
{
    pt::ptree& speciesA_db = calphad_db.get_child("SpeciesA");
    auto phase_db          = speciesA_db.get_child("PhaseL");
    g_species_phaseL_[0].initialize("L0", phase_db);
    readSublatticeStoichiometry(phase_db, sublattice_stoichiometry_phaseL_);

    phase_db = speciesA_db.get_child("PhaseA");
    g_species_phaseA_[0].initialize("A0", phase_db);
    readSublatticeStoichiometry(phase_db, sublattice_stoichiometry_phaseA_);

    pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
    g_species_phaseL_[1].initialize("L1", speciesB_db.get_child("PhaseL"));
    g_species_phaseA_[1].initialize("A1", speciesB_db.get_child("PhaseA"));

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

double CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    const CalphadDataType l0 = lmixPhase(0, pi, temperature);
    const CalphadDataType l1 = lmixPhase(1, pi, temperature);
    const CalphadDataType l2 = lmixPhase(2, pi, temperature);
    const CalphadDataType l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;
    CalphadDataType p;
    CalphadDataType q;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &g_species_phaseL_[0];
            p         = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseL_[0]);
            q = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseL_[1]);
            break;
        case PhaseIndex::phaseA:
            g_species = &g_species_phaseA_[0];
            p         = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseA_[0]);
            q = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseA_[1]);
            break;
        default:
            //            std::cerr <<
            //            "CALPHADFreeEnergyFunctionsBinary2Ph1Sl::"
            //                         "computeFreeEnergy(), undefined phase"
            //                      << "!!!" << std::endl;
            // abort();
            return 0.;
    }

    CalphadDataType ypp_A = (p + q) * conc[0] - p;

    double fe = (ypp_A * g_species[0].fenergy(temperature)
                    + (1. - ypp_A) * g_species[1].fenergy(temperature)
                    + CALPHADcomputeFMixBinary(l0, l1, l2, l3, ypp_A)
                    + CALPHADcomputeFIdealMixBinary(
                          GASCONSTANT_R_JPKPMOL * temperature, ypp_A))
                / (p + q);

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

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    const CalphadDataType l0 = lmixPhase(0, pi, temperature);
    const CalphadDataType l1 = lmixPhase(1, pi, temperature);
    const CalphadDataType l2 = lmixPhase(2, pi, temperature);
    const CalphadDataType l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;
    CalphadDataType p;
    CalphadDataType q;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &g_species_phaseL_[0];
            p         = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseL_[0]);
            q = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseL_[1]);
            break;
        case PhaseIndex::phaseA:
            g_species = &g_species_phaseA_[0];
            p         = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseA_[0]);
            q = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseA_[1]);
            break;
        default:
            //            std::cerr <<
            //            "CALPHADFreeEnergyFunctionsBinary2Ph1Sl::"
            //                         "computeFreeEnergy(), undefined phase!!!"
            //                      << std::endl;
            // abort();
            return;
    }

    CalphadDataType ypp_A = (p + q) * conc[0] - p;

    double mu = (g_species[0].fenergy(temperature)
                    - g_species[1].fenergy(temperature))
                + CALPHADcomputeFMix_derivBinary(l0, l1, l2, l3, ypp_A)
                + CALPHADcomputeFIdealMix_derivBinary(
                      GASCONSTANT_R_JPKPMOL * temperature, ypp_A);

    deriv[0] = mu;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    double* d2fdc2)
{
    // assert(conc[0] >= 0.);
    // assert(conc[0] <= 1.);

    const CalphadDataType l0 = lmixPhase(0, pi, temp);
    const CalphadDataType l1 = lmixPhase(1, pi, temp);
    const CalphadDataType l2 = lmixPhase(2, pi, temp);
    const CalphadDataType l3 = lmixPhase(3, pi, temp);
    const double rt          = GASCONSTANT_R_JPKPMOL * temp;

    CalphadDataType p;
    CalphadDataType q;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            p = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseL_[0]);
            q = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseL_[1]);
            break;
        case PhaseIndex::phaseA:
            p = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseA_[0]);
            q = static_cast<CalphadDataType>(
                sublattice_stoichiometry_phaseA_[1]);
            break;
        default:
            return;
    }

    CalphadDataType ypp_A = (p + q) * conc[0] - p;

    d2fdc2[0] = (CALPHADcomputeFMix_deriv2Binary(l0, l1, l2, l3, ypp_A)
                 + CALPHADcomputeFIdealMix_deriv2Binary(rt, ypp_A));
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computeTdependentParameters(
    const double temperature, CalphadDataType* Lmix_L, CalphadDataType* Lmix_A,
    CalphadDataType* fA, CalphadDataType* fB)
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
bool CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
    std::cerr << "CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computeCeqT() not "
                 "implemented"
              << std::endl;
    return false;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc,
    double& fl, double& fa)
{
    // std::cout<<"CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computePhasesFreeEnergies()"<<endl;

    double c[2] = { conc, conc };

    // evaluate temperature dependent parameters
    CalphadDataType fA[2];
    CalphadDataType fB[2];

    CalphadDataType Lmix_L[4];
    CalphadDataType Lmix_A[4];
    computeTdependentParameters(temperature, Lmix_L, Lmix_A, fA, fB);

    double RTinv = 1.0 / (GASCONSTANT_R_JPKPMOL * temperature);

    // Get the sublattice stoichiometry
    int p[2];
    int q[2];
    p[0] = sublattice_stoichiometry_phaseL_[0];
    p[1] = sublattice_stoichiometry_phaseA_[0];
    q[0] = sublattice_stoichiometry_phaseL_[1];
    q[1] = sublattice_stoichiometry_phaseA_[1];

    CALPHADConcSolverBinary2Ph1Sl solver;
    solver.setup(conc, hphi[0], RTinv, Lmix_L, Lmix_A, fA, fB, p, q);
    int ret = solver.ComputeConcentration(
        c, newton_tol_, newton_maxits_, newton_alpha_);
    if (ret < 0)
    {
#if 0
        std::cerr << "ERROR in "
                     "CALPHADFreeEnergyFunctionsBinary2Ph1Sl::"
                     "computePhasesFreeEnergies()"
                     " ---"
                  << "conc=" << conc << ", hphi=" << hphi << std::endl;
        abort();
#endif
    }

    // assert(c[0] >= 0.);
    fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);

    // assert(c[1] >= 0.);
    fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------

int CALPHADFreeEnergyFunctionsBinary2Ph1Sl::computePhaseConcentrations(
    const double temperature, const double* const conc, const double* const phi,
    double* x)
{
    // assert(x[0] >= 0.);
    // assert(x[1] >= 0.);
    // assert(x[0] <= 1.);
    // assert(x[1] <= 1.);

    const double RTinv = 1.0 / (GASCONSTANT_R_JPKPMOL * temperature);

    CalphadDataType fA[2];
    CalphadDataType fB[2];
    CalphadDataType Lmix_L[4];
    CalphadDataType Lmix_A[4];

    computeTdependentParameters(temperature, Lmix_L, Lmix_A, fA, fB);

    const double hphi = interp_func(conc_interp_func_type_, phi[0]);

    // Get the sublattice stoichiometry
    int p[2];
    int q[2];
    p[0] = sublattice_stoichiometry_phaseL_[0];
    p[1] = sublattice_stoichiometry_phaseA_[0];
    q[0] = sublattice_stoichiometry_phaseL_[1];
    q[1] = sublattice_stoichiometry_phaseA_[1];

    // conc could be outside of [0.,1.] in a trial step
    double c0 = conc[0] >= 0. ? conc[0] : 0.;
    c0        = c0 <= 1. ? c0 : 1.;
    // solve system of equations to find (cl,cs) given c0 and hphi
    // x: initial guess and solution
    CALPHADConcSolverBinary2Ph1Sl solver;
    solver.setup(c0, hphi, RTinv, Lmix_L, Lmix_A, fA, fB, p, q);
    int ret = solver.ComputeConcentration(
        x, newton_tol_, newton_maxits_, newton_alpha_);
#if 0
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "CALPHADFreeEnergyFunctionsBinary2Ph1Sl::"
                     "computePhaseConcentrations() "
                     "failed for conc="
                  << conc[0] << ", hphi=" << hphi[0] << std::endl;
        abort();
    }
#endif

    return ret;
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

//-----------------------------------------------------------------------

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::energyVsPhiAndC(
    const double temperature, const double* const ceq, const bool found_ceq,
    const double phi_well_scale, const int npts_phi, const int npts_c)
{
    std::cout << "CALPHADFreeEnergyFunctionsBinary2Ph1Sl::energyVsPhiAndC()..."
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
    std::cout << "CALPHADFreeEnergyFunctionsBinary2Ph1Sl: Use slope: " << slopec
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
void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::printEnergyVsPhiHeader(
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

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::printEnergyVsPhi(
    const double* const conc, const double temperature,
    const double phi_well_scale, const int npts, const double slopec,
    std::ostream& os)
{
    // std::cout <<
    // "CALPHADFreeEnergyFunctionsBinary2Ph1Sl::printEnergyVsPhi()..." <<
    // std::endl;
    const double dphi = 1.0 / (double)(npts - 1);

    // os << "# phi     f(phi)     for c=" << conc
    //           << " and T=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double phi = i * dphi;

        double e       = fchem(&phi, conc, temperature);
        const double w = phi_well_scale * well_func(phi);

        os << e + w + slopec * conc[0] << std::endl;
    }
    // os << std::endl;
}

//=======================================================================
// compute free energy in [J/mol]
double CALPHADFreeEnergyFunctionsBinary2Ph1Sl::fchem(
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

    double e = (1.0 - hfphi) * fl + hfphi * fa;

    return e;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::printEnergyVsComposition(
    const double temperature, std::ostream& os, const double cmin,
    const double cmax, const int npts)
{
    const double dc = (cmax - cmin) / (double)(npts - 1);

    os << "#phi=0" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        const double phi = 0.0;

        double e = fchem(&phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc + cmin;

        const double phi = 1.0;

        double e = fchem(&phi, &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary2Ph1Sl::preRunDiagnostics(
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
