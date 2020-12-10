#include "CALPHADFreeEnergyFunctionsBinary.h"
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
    const ConcInterpolationType conc_interp_func_type,
    const bool with_third_phase)
    : ceq_l_(-1.),
      ceq_a_(-1.),
      ceq_b_(-1.),
      energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type),
      with_third_phase_(with_third_phase),
      fenergy_diag_filename_("energy.vtk")
{
    readParameters(calphad_db);

    solver_ = new CALPHADConcentrationSolverBinary(with_third_phase_);

    if (newton_db) readNewtonparameters(newton_db.get());
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::readNewtonparameters(
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

void CALPHADFreeEnergyFunctionsBinary::readParameters(pt::ptree& calphad_db)
{
    pt::ptree& species0_db = calphad_db.get_child("SpeciesA");
    std::string dbnameL("PhaseL");
    g_species_phaseL_[0].initialize("L0", species0_db.get_child(dbnameL));
    std::string dbnameA("PhaseA");
    g_species_phaseA_[0].initialize("A0", species0_db.get_child(dbnameA));
    std::string dbnameB("PhaseB");
    if (with_third_phase_)
    {
        g_species_phaseB_[0].initialize("B0", species0_db.get_child(dbnameB));
    }

    pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
    g_species_phaseL_[1].initialize("L1", speciesB_db.get_child(dbnameL));
    g_species_phaseA_[1].initialize("A1", speciesB_db.get_child(dbnameA));
    if (with_third_phase_)
    {
        g_species_phaseB_[1].initialize("B1", speciesB_db.get_child(dbnameB));
    }

    // read Lmix coefficients
    std::string dbnamemixL("LmixPhaseL");
    pt::ptree Lmix0_db = calphad_db.get_child(dbnamemixL);
    readLmixBinary(Lmix0_db, LmixPhaseL_);

    std::string dbnamemixA("LmixPhaseA");
    pt::ptree Lmix1_db = calphad_db.get_child(dbnamemixA);
    readLmixBinary(Lmix1_db, LmixPhaseA_);

    if (with_third_phase_)
    {
        std::string dbnamemixB("LmixPhaseB");
        pt::ptree Lmix2_db = calphad_db.get_child(dbnamemixB);
        readLmixBinary(Lmix2_db, LmixPhaseB_);
    }

    // print database just read
    // std::clog << "CALPHAD database..." << std::endl;
    // pt::write_json(std::clog, calphad_db);
}

//-----------------------------------------------------------------------

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
        case PhaseIndex::phaseB:
            g_species = &g_species_phaseB_[0];
            break;
        default:
            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
                         "computeFreeEnergy(), undefined phase"
                      << "!!!" << std::endl;
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
        case PhaseIndex::phaseB:
            g_species = &g_species_phaseB_[0];
            break;
        default:
            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
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
    std::vector<double>& d2fdc2)
{
    assert(conc[0] >= 0.);
    assert(conc[0] <= 1.);

    const double l0 = lmixPhase(0, pi, temp);
    const double l1 = lmixPhase(1, pi, temp);
    const double l2 = lmixPhase(2, pi, temp);
    const double l3 = lmixPhase(3, pi, temp);
    const double rt = gas_constant_R_JpKpmol * temp;

    d2fdc2[0] = (CALPHADcomputeFMix_deriv2Binary(l0, l1, l2, l3, conc[0])
                 + CALPHADcomputeFIdealMix_deriv2Binary(rt, conc[0]));
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computeParametersForSolvers(
    const double temperature, double* Lmix_L, double* Lmix_A, double* Lmix_B,
    double* fA, double* fB, const PhaseIndex* pis, const int nphases)
{
    // loop over phases
    for (short i = 0; i < nphases; i++)
    {
        switch (pis[i])
        {
            case PhaseIndex::phaseL:
                fA[i]     = g_species_phaseL_[0].fenergy(temperature);
                fB[i]     = g_species_phaseL_[1].fenergy(temperature);
                Lmix_L[0] = lmixPhase(0, pis[i], temperature);
                Lmix_L[1] = lmixPhase(1, pis[i], temperature);
                Lmix_L[2] = lmixPhase(2, pis[i], temperature);
                Lmix_L[3] = lmixPhase(3, pis[i], temperature);
                break;

            case PhaseIndex::phaseA:
                fA[i]     = g_species_phaseA_[0].fenergy(temperature);
                fB[i]     = g_species_phaseA_[1].fenergy(temperature);
                Lmix_A[0] = lmixPhase(0, pis[i], temperature);
                Lmix_A[1] = lmixPhase(1, pis[i], temperature);
                Lmix_A[2] = lmixPhase(2, pis[i], temperature);
                Lmix_A[3] = lmixPhase(3, pis[i], temperature);
                break;

            case PhaseIndex::phaseB:
                fA[i]     = g_species_phaseB_[0].fenergy(temperature);
                fB[i]     = g_species_phaseB_[1].fenergy(temperature);
                Lmix_B[0] = lmixPhase(0, pis[i], temperature);
                Lmix_B[1] = lmixPhase(1, pis[i], temperature);
                Lmix_B[2] = lmixPhase(2, pis[i], temperature);
                Lmix_B[3] = lmixPhase(3, pis[i], temperature);
                break;

            default:
                std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
                             "computeParametersForSolvers: Undefined phase"
                          << std::endl;
        }
    }
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyFunctionsBinary::computeCeqT(const double temperature,
    const PhaseIndex pi0, const PhaseIndex pi1, double* ceq, const int maxits,
    const bool verbose)
{
    if (verbose)
        std::cout << "CALPHADFreeEnergyFunctionsBinary::computeCeqT()"
                  << std::endl;
    assert(temperature > 0.);

    // evaluate temperature dependent parameters
    double fA[3];
    double fB[3];

    double Lmix_L[4];
    double Lmix_A[4];
    double Lmix_B[4];

    PhaseIndex pis[2] = { pi0, pi1 };
    computeParametersForSolvers(
        temperature, Lmix_L, Lmix_A, Lmix_B, fA, fB, pis, 2);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    CALPHADEqConcSolverBinary eq_solver;
    eq_solver.SetMaxIterations(maxits);

    int ret = eq_solver.ComputeConcentration(
        ceq, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB);

    if (ret >= 0)
    {
        if (verbose)
        {
            std::cout << "CALPHAD, c_eq phase0=" << ceq[0] << std::endl;
            std::cout << "CALPHAD, c_eq phase1=" << ceq[1] << std::endl;
        }

        if (pi1 == PhaseIndex::phaseB)
        {
            ceq_b_ = ceq[1];
            ceq_l_ = 0.5 * (ceq[0] + ceq_l_);
        }
        else
        {
            ceq_l_ = ceq[0];
            ceq_a_ = ceq[1];
        }
    }
    else
    {
        std::cout << "CALPHADFreeEnergyFunctionsBinary, WARNING: ceq "
                     "computation did not converge"
                  << std::endl;
    }

    return (ret >= 0);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies(
    const double temperature, const double hphi, const double heta,
    const double conc, double& fl, double& fa, double& fb)
{
    // std::cout<<"CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;

    const int N = with_third_phase_ ? 3 : 2;

    double* c = new double[N];
    for (int ii = 0; ii < N; ii++)
    {
        c[ii] = conc;
    }
    // std::cout<<"ceq_l_="<<ceq_l_<<endl;
    // std::cout<<"ceq_a_="<<ceq_a_<<endl;
    // std::cout<<"ceq_b_="<<ceq_b_<<endl;
    if (ceq_l_ >= 0.) c[0] = ceq_l_;
    if (ceq_a_ >= 0.) c[1] = ceq_a_;
    if (ceq_b_ >= 0.) c[2] = ceq_b_;

    // evaluate temperature dependent parameters
    double fA[3];
    double fB[3];

    double Lmix_L[4];
    double Lmix_A[4];
    double Lmix_B[4];
    const PhaseIndex pis[3]
        = { PhaseIndex::phaseL, PhaseIndex::phaseA, PhaseIndex::phaseB };
    computeParametersForSolvers(
        temperature, Lmix_L, Lmix_A, Lmix_B, fA, fB, pis, 3);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);
    int ret      = solver_->ComputeConcentration(
        c, conc, hphi, heta, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB);

    if (ret < 0)
    {
        if (with_third_phase_)
        {
            std::cerr << "ceq_l_=" << ceq_l_ << std::endl;
            std::cerr << "ceq_a_=" << ceq_a_ << std::endl;
            std::cerr << "ceq_b_=" << ceq_b_ << std::endl;
            std::cerr << "ERROR in "
                         "CALPHADFreeEnergyFunctionsBinary::"
                         "computePhasesFreeEnergies()"
                         " ---"
                      << "conc=" << conc << ", hphi=" << hphi
                      << ", heta=" << heta << std::endl;
        }
        else
        {
            std::cerr << "ERROR in "
                         "CALPHADFreeEnergyFunctionsBinary::"
                         "computePhasesFreeEnergies()"
                         " ---"
                      << "conc=" << conc << ", hphi=" << hphi << std::endl;
        }
        abort();
    }

    assert(c[0] >= 0.);
    fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);

    assert(c[1] >= 0.);
    fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);

    fb = 0.;
    if (with_third_phase_)
    {
        assert(c[2] >= 0.);
        assert(c[2] <= 1.);
        fb = computeFreeEnergy(temperature, &c[2], PhaseIndex::phaseB, false);
    }
    delete[] c;
}

//-----------------------------------------------------------------------

int CALPHADFreeEnergyFunctionsBinary::computePhaseConcentrations(
    const double temperature, const double* const conc, const double phi,
    const double eta, double* x)
{
    assert(x[0] >= 0.);
    assert(x[1] >= 0.);
    assert(x[0] <= 1.);
    assert(x[1] <= 1.);

    const double conc0 = conc[0];

    const double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    double fA[3];
    double fB[3];
    double Lmix_L[4];
    double Lmix_A[4];
    double Lmix_B[4];

    const int nphases = with_third_phase_ ? 3 : 2;
    PhaseIndex* pis   = new PhaseIndex[nphases];
    pis[0]            = PhaseIndex::phaseL;
    pis[1]            = PhaseIndex::phaseA;
    if (with_third_phase_) pis[2] = PhaseIndex::phaseB;
    computeParametersForSolvers(
        temperature, Lmix_L, Lmix_A, Lmix_B, fA, fB, pis, nphases);

    const char interp_func_type = concInterpChar(conc_interp_func_type_);
    const double hphi           = interp_func(phi, interp_func_type);

    double heta = 0.0;
    // std::cout<<"ceq_a_="<<ceq_a_<<endl;
    // x[0] = ( ceq_l_>=0. ) ? ceq_l_ : 0.5;
    // x[1] = ( ceq_a_>=0. ) ? ceq_a_ : 0.5;

    if (with_third_phase_)
    {
        assert(x[2] >= 0.);
        assert(x[2] <= 1.);
        // x[2] = ( ceq_b_>=0. ) ? ceq_b_ : 0.5;

        heta = interp_func(eta, interp_func_type);
    }

    // conc could be outside of [0.,1.] in a trial step
    double c0 = conc[0] >= 0. ? conc[0] : 0.;
    c0        = c0 <= 1. ? c0 : 1.;
    int ret   = solver_->ComputeConcentration(
        x, c0, hphi, heta, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB);
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "CALPHADFreeEnergyFunctionsBinary::"
                     "computePhaseConcentrations() "
                     "failed for conc="
                  << conc0 << ", hphi=" << hphi << ", heta=" << heta
                  << std::endl;
        abort();
    }

    return ret;
}

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
    const double eta  = 0.0;

    // os << "# phi     f(phi)     for c=" << conc
    //           << " eta=" << eta
    //           << " and T=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double phi = i * dphi;

        double e       = fchem(phi, eta, conc, temperature);
        const double w = phi_well_scale * well_func(phi);

        os << e + w + slopec * conc[0] << std::endl;
    }
    // os << std::endl;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsEta(
    const double* const conc, const double temperature,
    const double eta_well_scale, const int npts, const double slopec,
    std::ostream& os)
{
    // std::cout << "CALPHADFreeEnergyFunctionsBinary::printEnergyVsEta()..." <<
    // std::endl;
    const double deta = 1.0 / (double)(npts - 1);
    const double phi  = 1.0;

    for (int i = 0; i < npts; i++)
    {
        const double eta = i * deta;

        double e = fchem(phi, eta, conc, temperature);

        double w = eta_well_scale * well_func(eta);

        os << e + w + slopec * conc[0] << std::endl;
    }
    // os << std::endl;
}

//=======================================================================
// compute free energy in [J/mol]
double CALPHADFreeEnergyFunctionsBinary::fchem(const double phi,
    const double eta, const double* const conc, const double temperature)
{
    const char interp_func_type = concInterpChar(conc_interp_func_type_);
    const double hcphi          = interp_func(phi, interp_func_type);
    double heta                 = 0.0;
    if (with_third_phase_)
    {
        heta = interp_func(eta, interp_func_type);
    }

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;
    double fb        = 0.;
    if ((phi > tol) & (phi < (1. - tol)))
    {
        computePhasesFreeEnergies(
            temperature, hcphi, heta, conc[0], fl, fa, fb);
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
            if (with_third_phase_)
                fb = computeFreeEnergy(temperature, conc, PhaseIndex::phaseB);
        }
    }

    const char interpf = energyInterpChar(energy_interp_func_type_);
    const double hfphi = interp_func(phi, interpf);
    double e = (1.0 - hfphi) * fl + hfphi * ((1.0 - heta) * fa + heta * fb);

    return e;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsComposition(
    const double temperature, const int npts)
{
    std::ofstream os("FvsC.dat", std::ios::out);

    const double dc = 1.0 / (double)(npts - 1);

    os << "#phi=0" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        double e = fchem(0., 0., &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        double e = fchem(1., 0., &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }

    if (with_third_phase_)
    {
        os << std::endl;

        os << "#eta=1" << std::endl;
        for (int i = 0; i < npts; i++)
        {
            const double conc = i * dc;

            double e = fchem(1., 1., &conc, temperature);
            os << conc << "\t" << e << std::endl;
        }
    }
}
}
