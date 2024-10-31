#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADConcSolverBinaryThreePhase.h"
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

CALPHADFreeEnergyFunctionsBinaryThreePhase::
    CALPHADFreeEnergyFunctionsBinaryThreePhase(pt::ptree& calphad_db,
        boost::optional<pt::ptree&> newton_db,
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
void CALPHADFreeEnergyFunctionsBinaryThreePhase::readNewtonparameters(
    pt::ptree& newton_db)
{
    newton_tol_     = newton_db.get<double>("tol", newton_tol_);
    newton_alpha_   = newton_db.get<double>("alpha", newton_alpha_);
    newton_maxits_  = newton_db.get<int>("max_its", newton_maxits_);
    newton_verbose_ = newton_db.get<bool>("verbose", newton_verbose_);

    std::cout << "Thermo4PFM Newton Solver Parameters:" << std::endl;
    std::cout << "Tolerance: " << newton_tol_ << std::endl;
    std::cout << "Alpha: " << newton_alpha_ << std::endl;
    std::cout << "Max interations: " << newton_maxits_ << std::endl;
    std::cout << "Verbose: " << newton_verbose_ << std::endl;
}

//=======================================================================
void CALPHADFreeEnergyFunctionsBinaryThreePhase::readParameters(
    pt::ptree& calphad_db)
{
    // print database to be read
    // std::clog << "CALPHAD database..." << std::endl;
    // pt::write_json(std::clog, calphad_db);

    pt::ptree& speciesA_db = calphad_db.get_child("SpeciesA");
    std::string dbnameL("PhaseL");
    g_species_phaseL_[0].initialize("L0", speciesA_db.get_child(dbnameL));
    // assert(checkSingleSublattice(speciesA_db.get_child(dbnameL)));
    std::string dbnameA("PhaseA");
    g_species_phaseA_[0].initialize("A0", speciesA_db.get_child(dbnameA));
    // assert(checkSingleSublattice(speciesA_db.get_child(dbnameA)));
    std::string dbnameB("PhaseB");
    g_species_phaseB_[0].initialize("B0", speciesA_db.get_child(dbnameB));
    // assert(checkSingleSublattice(speciesA_db.get_child(dbnameB)));

    pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
    g_species_phaseL_[1].initialize("L1", speciesB_db.get_child(dbnameL));
    g_species_phaseA_[1].initialize("A1", speciesB_db.get_child(dbnameA));
    g_species_phaseB_[1].initialize("B1", speciesB_db.get_child(dbnameB));

    // read Lmix coefficients
    std::string dbnamemixL("LmixPhaseL");
    pt::ptree Lmix0_db = calphad_db.get_child(dbnamemixL);
    readLmixBinary(Lmix0_db, LmixPhaseL_);

    std::string dbnamemixA("LmixPhaseA");
    pt::ptree Lmix1_db = calphad_db.get_child(dbnamemixA);
    readLmixBinary(Lmix1_db, LmixPhaseA_);

    std::string dbnamemixB("LmixPhaseB");
    pt::ptree Lmix2_db = calphad_db.get_child(dbnamemixB);
    readLmixBinary(Lmix2_db, LmixPhaseB_);
}

//-----------------------------------------------------------------------
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double CALPHADFreeEnergyFunctionsBinaryThreePhase::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    const CalphadDataType l0 = lmixPhase(0, pi, temperature);
    const CalphadDataType l1 = lmixPhase(1, pi, temperature);
    const CalphadDataType l2 = lmixPhase(2, pi, temperature);
    const CalphadDataType l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &g_species_phaseL_[0];
            break;
        case PhaseIndex::phaseA:
            g_species = &g_species_phaseA_[0];
            break;
        case PhaseIndex::phaseB:
            g_species = &g_species_phaseB_[0];
            break;
        default:
            //            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
            //                         "computeFreeEnergy(), undefined phase"
            //                      << "!!!" << std::endl;
            // abort();
            return 0.;
    }

    double fe = conc[0] * g_species[0].fenergy(temperature)
                + (1. - conc[0]) * g_species[1].fenergy(temperature)
                + CALPHADcomputeFMixBinary(l0, l1, l2, l3, conc[0])
                + CALPHADcomputeFIdealMixBinary(
                      GASCONSTANT_R_JPKPMOL * temperature, conc[0]);

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

void CALPHADFreeEnergyFunctionsBinaryThreePhase::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    const CalphadDataType l0 = lmixPhase(0, pi, temperature);
    const CalphadDataType l1 = lmixPhase(1, pi, temperature);
    const CalphadDataType l2 = lmixPhase(2, pi, temperature);
    const CalphadDataType l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &g_species_phaseL_[0];
            break;
        case PhaseIndex::phaseA:
            g_species = &g_species_phaseA_[0];
            break;
        case PhaseIndex::phaseB:
            g_species = &g_species_phaseB_[0];
            break;
        default:
            //            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
            //                         "computeFreeEnergy(), undefined phase!!!"
            //                      << std::endl;
            // abort();
            return;
    }

    double mu = (g_species[0].fenergy(temperature)
                    - g_species[1].fenergy(temperature))
                + CALPHADcomputeFMix_derivBinary(l0, l1, l2, l3, conc[0])
                + CALPHADcomputeFIdealMix_derivBinary(
                      GASCONSTANT_R_JPKPMOL * temperature, conc[0]);

    deriv[0] = mu;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinaryThreePhase::
    computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2)
{
    const CalphadDataType l0 = lmixPhase(0, pi, temp);
    const CalphadDataType l1 = lmixPhase(1, pi, temp);
    const CalphadDataType l2 = lmixPhase(2, pi, temp);
    const CalphadDataType l3 = lmixPhase(3, pi, temp);
    const double rt          = GASCONSTANT_R_JPKPMOL * temp;

    d2fdc2[0] = (CALPHADcomputeFMix_deriv2Binary(l0, l1, l2, l3, conc[0])
                 + CALPHADcomputeFIdealMix_deriv2Binary(rt, conc[0]));
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinaryThreePhase::computeTdependentParameters(
    const double temperature, CalphadDataType* Lmix_L, CalphadDataType* Lmix_A,
    CalphadDataType* Lmix_B, CalphadDataType* fA, CalphadDataType* fB)
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

    fA[2]     = g_species_phaseB_[0].fenergy(temperature);
    fB[2]     = g_species_phaseB_[1].fenergy(temperature);
    Lmix_B[0] = lmixPhase(0, PhaseIndex::phaseB, temperature);
    Lmix_B[1] = lmixPhase(1, PhaseIndex::phaseB, temperature);
    Lmix_B[2] = lmixPhase(2, PhaseIndex::phaseB, temperature);
    Lmix_B[3] = lmixPhase(3, PhaseIndex::phaseB, temperature);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinaryThreePhase::computePhasesFreeEnergies(
    const double temperature, const double* const hphi, const double conc,
    double& fl, double& fa, double& fb)
{
    // std::cout<<"CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;

    double c[3] = { conc, conc, conc };

    int ret = computePhaseConcentrations(temperature, &conc, hphi, c);

    if (ret < 0)
    {
#if 0
        std::cerr << "ERROR in "
                     "CALPHADFreeEnergyFunctionsBinaryThreePhase::"
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

int CALPHADFreeEnergyFunctionsBinaryThreePhase::computePhaseConcentrations(
    const double temperature, const double* const conc, const double* const phi,
    double* x)
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

    // solve system of equations to find (cl,cs) given conc[0] and hphi
    // x: initial guess and solution
    CALPHADConcSolverBinaryThreePhase solver;
    solver.setup(
        conc[0], hphi0, hphi1, hphi2, RT, Lmix_L, Lmix_A, Lmix_B, fA, fB);
    int ret = solver.ComputeConcentration(
        x, newton_tol_, newton_maxits_, newton_alpha_);
#if 0
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "CALPHADFreeEnergyFunctionsBinary::"
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

double CALPHADFreeEnergyFunctionsBinaryThreePhase::fchem(
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

    double e = hfphi[0] * fl + hfphi[1] * fa + hfphi[2] * fb;

    return e;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinaryThreePhase::printEnergyVsComposition(
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

void CALPHADFreeEnergyFunctionsBinaryThreePhase::preRunDiagnostics(
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
