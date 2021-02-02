#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "CALPHADFunctions.h"
#include "PhysicalConstants.h"
#include "functions.h"
#include "well_functions.h"

#include <boost/property_tree/json_parser.hpp>

#include <cmath>
#include <iomanip>
#include <string>

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

void readLmixTernaryParameters(pt::ptree& Lmix_db, double LmixABC[3][2])
{
    // L0
    {
        auto child = Lmix_db.get_child_optional("L0");
        if (child)
        {
            int i = 0;
            for (pt::ptree::value_type& v : Lmix_db.get_child("L0"))
            {
                LmixABC[0][i] = v.second.get_value<double>();
                i++;
            }
        }
    }
    // L1
    {
        auto child = Lmix_db.get_child_optional("L1");
        if (child)
        {
            int i = 0;
            for (pt::ptree::value_type& v : Lmix_db.get_child("L1"))
            {
                LmixABC[1][i] = v.second.get_value<double>();
                i++;
            }
        }
    }
    // L2
    {
        auto child = Lmix_db.get_child_optional("L2");
        if (child)
        {
            int i = 0;
            for (pt::ptree::value_type& v : Lmix_db.get_child("L2"))
            {
                LmixABC[2][i] = v.second.get_value<double>();
                i++;
            }
        }
    }
}
CALPHADFreeEnergyFunctionsTernary::CALPHADFreeEnergyFunctionsTernary(
    pt::ptree& calphad_db, boost::optional<pt::ptree&> newton_db,
    const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type)
    : energy_interp_func_type_(energy_interp_func_type),
      conc_interp_func_type_(conc_interp_func_type),
      fenergy_diag_filename_("energy.vtk"),
      newton_tol_(1.e-8),
      newton_alpha_(1.),
      newton_maxits_(20),
      newton_verbose_(false)

{
    readParameters(calphad_db);

    if (newton_db) readNewtonparameters(newton_db.get());
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::readNewtonparameters(
    pt::ptree& newton_db)
{
    newton_tol_     = newton_db.get<double>("tol", newton_tol_);
    newton_alpha_   = newton_db.get<double>("alpha", newton_alpha_);
    newton_maxits_  = newton_db.get<int>("max_its", newton_maxits_);
    newton_verbose_ = newton_db.get<bool>("verbose", newton_verbose_);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::readParameters(pt::ptree& calphad_db)
{
    pt::ptree& species0_db = calphad_db.get_child("SpeciesA");
    std::string name       = species0_db.get<std::string>("name", "unknown");
    std::string dbnameL("PhaseL");
    g_species_phaseL_[0].initialize(name, species0_db.get_child(dbnameL));
    std::string dbnameA("PhaseA");
    g_species_phaseA_[0].initialize(name, species0_db.get_child(dbnameA));

    pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
    name                   = speciesB_db.get<std::string>("name", "unknown");
    g_species_phaseL_[1].initialize(name, speciesB_db.get_child(dbnameL));
    g_species_phaseA_[1].initialize(name, speciesB_db.get_child(dbnameA));

    pt::ptree& speciesC_db = calphad_db.get_child("SpeciesC");
    name                   = speciesC_db.get<std::string>("name", "unknown");
    g_species_phaseL_[2].initialize(name, speciesC_db.get_child(dbnameL));
    g_species_phaseA_[2].initialize(name, speciesC_db.get_child(dbnameA));

    // read Lmix coefficients

    // AB
    {
        pt::ptree& Lmix0_db = calphad_db.get_child("LmixABPhaseL");
        readLmixBinary(Lmix0_db, LmixABPhaseL_);

        pt::ptree& Lmix1_db = calphad_db.get_child("LmixABPhaseA");
        readLmixBinary(Lmix1_db, LmixABPhaseA_);
    }

    // AC
    {
        pt::ptree& Lmix0_db = calphad_db.get_child("LmixACPhaseL");
        readLmixBinary(Lmix0_db, LmixACPhaseL_);

        pt::ptree& Lmix1_db = calphad_db.get_child("LmixACPhaseA");
        readLmixBinary(Lmix1_db, LmixACPhaseA_);
    }

    // BC
    {
        pt::ptree& Lmix0_db = calphad_db.get_child("LmixBCPhaseL");
        readLmixBinary(Lmix0_db, LmixBCPhaseL_);

        pt::ptree& Lmix1_db = calphad_db.get_child("LmixBCPhaseA");
        readLmixBinary(Lmix1_db, LmixBCPhaseA_);
    }

    // ABC
    {
        // default values
        LmixABCPhaseL_[0][0] = 0.0;
        LmixABCPhaseL_[0][1] = 0.0;
        LmixABCPhaseL_[1][0] = 0.0;
        LmixABCPhaseL_[1][1] = 0.0;
        LmixABCPhaseL_[2][0] = 0.0;
        LmixABCPhaseL_[2][1] = 0.0;

        std::string dbnamemixL("LmixABCPhaseL");
        if (calphad_db.get_child_optional(dbnamemixL))
        {
            pt::ptree& Lmix0_db = calphad_db.get_child(dbnamemixL);
            readLmixTernaryParameters(Lmix0_db, LmixABCPhaseL_);
        }

        assert(LmixABCPhaseL_[0][0] == LmixABCPhaseL_[0][0]);
        assert(LmixABCPhaseL_[0][1] == LmixABCPhaseL_[0][1]);
        assert(LmixABCPhaseL_[1][0] == LmixABCPhaseL_[1][0]);
        assert(LmixABCPhaseL_[1][1] == LmixABCPhaseL_[1][1]);
        assert(LmixABCPhaseL_[2][0] == LmixABCPhaseL_[2][0]);
        assert(LmixABCPhaseL_[2][1] == LmixABCPhaseL_[2][1]);

        std::string dbnamemixA("LmixABCPhaseA");
        // default values
        LmixABCPhaseA_[0][0] = 0.0;
        LmixABCPhaseA_[0][1] = 0.0;
        LmixABCPhaseA_[1][0] = 0.0;
        LmixABCPhaseA_[1][1] = 0.0;
        LmixABCPhaseA_[2][0] = 0.0;
        LmixABCPhaseA_[2][1] = 0.0;
        if (calphad_db.get_child_optional(dbnamemixA))
        {
            pt::ptree& Lmix1_db = calphad_db.get_child(dbnamemixA);
            readLmixTernaryParameters(Lmix1_db, LmixABCPhaseA_);
        }

        assert(LmixABCPhaseA_[0][0] == LmixABCPhaseA_[0][0]);
        assert(LmixABCPhaseA_[0][1] == LmixABCPhaseA_[0][1]);
        assert(LmixABCPhaseA_[1][0] == LmixABCPhaseA_[1][0]);
        assert(LmixABCPhaseA_[1][1] == LmixABCPhaseA_[1][1]);
        assert(LmixABCPhaseA_[2][0] == LmixABCPhaseA_[2][0]);
        assert(LmixABCPhaseA_[2][1] == LmixABCPhaseA_[2][1]);
    }

    // print database just read
    std::clog << "CALPHAD database..." << std::endl;
    // pt::write_json(std::clog, calphad_db);
}

//-----------------------------------------------------------------------

double CALPHADFreeEnergyFunctionsTernary::computeFreeEnergy(
    const double temperature, const double* conc, const PhaseIndex pi,
    const bool gp)
{
    const double conc0 = conc[0];
    const double conc1 = conc[1];

    double lAB[4]
        = { lmix0ABPhase(pi, temperature), lmix1ABPhase(pi, temperature),
              lmix2ABPhase(pi, temperature), lmix3ABPhase(pi, temperature) };
    double lAC[4]
        = { lmix0ACPhase(pi, temperature), lmix1ACPhase(pi, temperature),
              lmix2ACPhase(pi, temperature), lmix3ACPhase(pi, temperature) };
    double lBC[4]
        = { lmix0BCPhase(pi, temperature), lmix1BCPhase(pi, temperature),
              lmix2BCPhase(pi, temperature), lmix3BCPhase(pi, temperature) };

    double lABC[3] = { lmix0ABCPhase(pi, temperature),
        lmix1ABCPhase(pi, temperature), lmix2ABCPhase(pi, temperature) };

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
            std::cout << "CALPHADFreeEnergyFunctionsTernary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
            abort();
            return 0.;
    }

    double conc2 = 1. - conc0 - conc1;
    double fe    = conc0 * g_species[0].fenergy(temperature)
                + conc1 * g_species[1].fenergy(temperature)
                + conc2 * g_species[2].fenergy(temperature)
                + CALPHADcomputeFMixTernary(lAB, lAC, lBC, lABC, conc0, conc1)
                + CALPHADcomputeFIdealMixTernary(
                      gas_constant_R_JpKpmol * temperature, conc0, conc1);

    // subtract -mu*c to get grand potential
    if (gp)
    {
        double deriv[2];
        computeDerivFreeEnergy(temperature, conc, pi, deriv);
        fe -= deriv[0] * conc0;
        fe -= deriv[1] * conc1;
    }

    return fe;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    double lAB[4]
        = { lmix0ABPhase(pi, temperature), lmix1ABPhase(pi, temperature),
              lmix2ABPhase(pi, temperature), lmix3ABPhase(pi, temperature) };
    double lAC[4]
        = { lmix0ACPhase(pi, temperature), lmix1ACPhase(pi, temperature),
              lmix2ACPhase(pi, temperature), lmix3ACPhase(pi, temperature) };
    double lBC[4]
        = { lmix0BCPhase(pi, temperature), lmix1BCPhase(pi, temperature),
              lmix2BCPhase(pi, temperature), lmix3BCPhase(pi, temperature) };
    double lABC[3] = { lmix0ABCPhase(pi, temperature),
        lmix1ABCPhase(pi, temperature), lmix2ABCPhase(pi, temperature) };

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
            std::cout << "CALPHADFreeEnergyFunctionsTernary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
            abort();
            return;
    }

    CALPHADcomputeFMix_derivTernary(
        lAB, lAC, lBC, lABC, conc[0], conc[1], deriv);

    deriv[0] += g_species[0].fenergy(temperature);
    deriv[0] -= g_species[2].fenergy(temperature);

    deriv[1] += g_species[1].fenergy(temperature);
    deriv[1] -= g_species[2].fenergy(temperature);

    double tmp[2];
    CALPHADcomputeFIdealMix_derivTernary(
        gas_constant_R_JpKpmol * temperature, conc[0], conc[1], tmp);
    deriv[0] += tmp[0];
    deriv[1] += tmp[1];
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    double* d2fdc2)
{
    assert(conc[0] >= 0.);
    assert(conc[0] <= 1.);
    assert(conc[1] >= 0.);
    assert(conc[1] <= 1.);

    double lAB[4]   = { lmix0ABPhase(pi, temp), lmix1ABPhase(pi, temp),
        lmix2ABPhase(pi, temp), lmix3ABPhase(pi, temp) };
    double lAC[4]   = { lmix0ACPhase(pi, temp), lmix1ACPhase(pi, temp),
        lmix2ACPhase(pi, temp), lmix3ACPhase(pi, temp) };
    double lBC[4]   = { lmix0BCPhase(pi, temp), lmix1BCPhase(pi, temp),
        lmix2BCPhase(pi, temp), lmix3BCPhase(pi, temp) };
    double lABC[3]  = { lmix0ABCPhase(pi, temp), lmix1ABCPhase(pi, temp),
        lmix2ABCPhase(pi, temp) };
    const double rt = gas_constant_R_JpKpmol * temp;

    double deriv1[4];
    CALPHADcomputeFIdealMix_deriv2Ternary(rt, conc[0], conc[1], &deriv1[0]);

    double deriv2[4];
    CALPHADcomputeFMix_deriv2Ternary(
        lAB, lAC, lBC, lABC, conc[0], conc[1], &deriv2[0]);

    d2fdc2[0] = deriv1[0] + deriv2[0];
    d2fdc2[1] = deriv1[1] + deriv2[1];
    d2fdc2[2] = deriv1[2] + deriv2[2];
    d2fdc2[3] = deriv1[3] + deriv2[3];
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::computeTdependentParameters(
    const double temperature, double* L_AB_L, double* L_AC_L, double* L_BC_L,
    double* L_ABC_L, double* L_AB_S, double* L_AC_S, double* L_BC_S,
    double* L_ABC_S, double* fA, double* fB, double* fC)
{
    fA[0] = g_species_phaseL_[0].fenergy(temperature);
    fB[0] = g_species_phaseL_[1].fenergy(temperature);
    fC[0] = g_species_phaseL_[2].fenergy(temperature);

    fA[1] = g_species_phaseA_[0].fenergy(temperature);
    fB[1] = g_species_phaseA_[1].fenergy(temperature);
    fC[1] = g_species_phaseA_[2].fenergy(temperature);

    L_AB_L[0]  = lmix0ABPhaseL(temperature);
    L_AB_L[1]  = lmix1ABPhaseL(temperature);
    L_AB_L[2]  = lmix2ABPhaseL(temperature);
    L_AB_L[3]  = lmix3ABPhaseL(temperature);
    L_AC_L[0]  = lmix0ACPhaseL(temperature);
    L_AC_L[1]  = lmix1ACPhaseL(temperature);
    L_AC_L[2]  = lmix2ACPhaseL(temperature);
    L_AC_L[3]  = lmix3ACPhaseL(temperature);
    L_BC_L[0]  = lmix0BCPhaseL(temperature);
    L_BC_L[1]  = lmix1BCPhaseL(temperature);
    L_BC_L[2]  = lmix2BCPhaseL(temperature);
    L_BC_L[3]  = lmix3BCPhaseL(temperature);
    L_ABC_L[0] = lmix0ABCPhaseL(temperature);
    L_ABC_L[1] = lmix1ABCPhaseL(temperature);
    L_ABC_L[2] = lmix2ABCPhaseL(temperature);

    L_AB_S[0]  = lmix0ABPhaseA(temperature);
    L_AB_S[1]  = lmix1ABPhaseA(temperature);
    L_AB_S[2]  = lmix2ABPhaseA(temperature);
    L_AB_S[3]  = lmix3ABPhaseA(temperature);
    L_AC_S[0]  = lmix0ACPhaseA(temperature);
    L_AC_S[1]  = lmix1ACPhaseA(temperature);
    L_AC_S[2]  = lmix2ACPhaseA(temperature);
    L_AC_S[3]  = lmix3ACPhaseA(temperature);
    L_BC_S[0]  = lmix0BCPhaseA(temperature);
    L_BC_S[1]  = lmix1BCPhaseA(temperature);
    L_BC_S[2]  = lmix2BCPhaseA(temperature);
    L_BC_S[3]  = lmix3BCPhaseA(temperature);
    L_ABC_S[0] = lmix0ABCPhaseA(temperature);
    L_ABC_S[1] = lmix1ABCPhaseA(temperature);
    L_ABC_S[2] = lmix2ABCPhaseA(temperature);

    assert(L_ABC_L[0] == L_ABC_L[0]);
    assert(fC[0] == fC[0]);
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyFunctionsTernary::computeCeqT(
    const double temperature, double* ceq, const int maxits, const bool verbose)
{
    if (verbose)
        std::cout << "CALPHADFreeEnergyFunctionsTernary::computeCeqT()"
                  << std::endl;
    assert(temperature > 0.);

    double L_AB_L[4];
    double L_AC_L[4];
    double L_BC_L[4];

    double L_ABC_L[3];

    double L_AB_S[4];
    double L_AC_S[4];
    double L_BC_S[4];

    double L_ABC_S[3];

    double fA[2];
    double fB[2];
    double fC[2];

    computeTdependentParameters(temperature, L_AB_L, L_AC_L, L_BC_L, L_ABC_L,
        L_AB_S, L_AC_S, L_BC_S, L_ABC_S, fA, fB, fC);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    CALPHADEqConcSolverTernary eq_solver;
    eq_solver.setup(RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S, L_AC_S, L_BC_S,
        L_ABC_L, L_ABC_S, fA, fB, fC);
    int ret = eq_solver.ComputeConcentration(ceq, newton_tol_, maxits);

    if (ret >= 0)
    {
        if (verbose)
        {
            std::cout << "CALPHAD, c0 phase0=" << ceq[0] << std::endl;
            std::cout << "CALPHAD, c1 phase0=" << ceq[1] << std::endl;
            std::cout << "CALPHAD, c0 phase1=" << ceq[2] << std::endl;
            std::cout << "CALPHAD, c1 phase1=" << ceq[3] << std::endl;
        }
    }
    else
    {
        std::cout << "CALPHADFreeEnergyFunctionsTernary, WARNING: ceq "
                     "computation did not converge"
                  << std::endl;
    }

    return (ret >= 0);
}

//=======================================================================

bool CALPHADFreeEnergyFunctionsTernary::computeCeqT(const double temperature,
    const double c0, const double c1, double* ceq, const int maxits,
    const bool verbose)
{
    assert(temperature > 0.);

    double L_AB_L[4];
    double L_AC_L[4];
    double L_BC_L[4];

    double L_ABC_L[3];

    double L_AB_S[4];
    double L_AC_S[4];
    double L_BC_S[4];

    double L_ABC_S[3];

    double fA[2];
    double fB[2];
    double fC[2];

    computeTdependentParameters(temperature, L_AB_L, L_AC_L, L_BC_L, L_ABC_L,
        L_AB_S, L_AC_S, L_BC_S, L_ABC_S, fA, fB, fC);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    CALPHADEqPhaseConcSolverTernary eq_solver;
    eq_solver.setup(c0, c1, RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S, L_AC_S,
        L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);
    int ret = eq_solver.ComputeConcentration(ceq, newton_tol_, newton_maxits_);
    if (ret >= 0)
    {
        if (verbose)
        {
            std::cout << "CALPHAD, c0 phase0=" << ceq[0] << std::endl;
            std::cout << "CALPHAD, c1 phase0=" << ceq[1] << std::endl;
            std::cout << "CALPHAD, c0 phase1=" << ceq[2] << std::endl;
            std::cout << "CALPHAD, c1 phase1=" << ceq[3] << std::endl;
        }
    }
    else
    {
        std::cout << "CALPHADFreeEnergyFunctionsTernary, WARNING: ceq "
                     "computation did not converge"
                  << std::endl;
    }

    return (ret >= 0);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::computePhasesFreeEnergies(
    const double temperature, const double hphi, const double conc0,
    const double conc1, double& fl, double& fa)
{
    // std::cout<<"CALPHADFreeEnergyFunctionsTernary::computePhasesFreeEnergies()"<<endl;

    double cauxilliary[4] = { conc0, conc1, conc0, conc1 };

    double L_AB_L[4];
    double L_AC_L[4];
    double L_BC_L[4];

    double L_ABC_L[3];

    double L_AB_S[4];
    double L_AC_S[4];
    double L_BC_S[4];

    double L_ABC_S[3];

    double fA[2];
    double fB[2];
    double fC[2];

    computeTdependentParameters(temperature, L_AB_L, L_AC_L, L_BC_L, L_ABC_L,
        L_AB_S, L_AC_S, L_BC_S, L_ABC_S, fA, fB, fC);

    assert(fC[0] == fC[0]);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);
    CALPHADConcSolverTernary solver;
    solver.setup(conc0, conc1, hphi, RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S,
        L_AC_S, L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);
    int ret
        = solver.ComputeConcentration(cauxilliary, newton_tol_, newton_maxits_);
    if (ret < 0)
    {
        std::cerr << "ERROR in "
                     "CALPHADFreeEnergyFunctionsTernary::"
                     "computePhasesFreeEnergies() "
                     "---"
                  << "conc0=" << conc0 << ", conc1=" << conc1
                  << ", hphi=" << hphi << std::endl;
        abort();
    }

    assert(conc0 >= 0.);
    double concl[2] = { cauxilliary[0], cauxilliary[1] };
    fl = computeFreeEnergy(temperature, &concl[0], PhaseIndex::phaseL, false);

    assert(conc1 >= 0.);
    double conca[2] = { cauxilliary[2], cauxilliary[3] };
    fa = computeFreeEnergy(temperature, &conca[0], PhaseIndex::phaseA, false);
}

//-----------------------------------------------------------------------
// output: x
int CALPHADFreeEnergyFunctionsTernary::computePhaseConcentrations(
    const double temperature, const double* const conc, const double phi,
    double* x)
{
    assert(conc[0] == conc[0]);
    assert(conc[1] == conc[1]);
    assert(x[0] >= 0.);
    assert(x[1] >= 0.);
    assert(x[0] <= 1.);
    assert(x[1] <= 1.);

    const double conc0 = conc[0];
    const double conc1 = conc[1];

    const double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    double L_AB_L[4];
    double L_AC_L[4];
    double L_BC_L[4];

    double L_ABC_L[3];

    double L_AB_S[4];
    double L_AC_S[4];
    double L_BC_S[4];

    double L_ABC_S[3];

    double fA[2];
    double fB[2];
    double fC[2];

    computeTdependentParameters(temperature, L_AB_L, L_AC_L, L_BC_L, L_ABC_L,
        L_AB_S, L_AC_S, L_BC_S, L_ABC_S, fA, fB, fC);
    assert(fC[0] == fC[0]);

    const double hphi = interp_func(conc_interp_func_type_, phi);

    // conc could be outside of [0.,1.] in a trial step
    double c0 = conc0 >= 0. ? conc0 : 0.;
    c0        = c0 <= 1. ? c0 : 1.;
    double c1 = conc1 >= 0. ? conc1 : 0.;
    c1        = c1 <= 1. ? c1 : 1.;

    CALPHADConcSolverTernary solver;
    solver.setup(c0, c1, hphi, RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S, L_AC_S,
        L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);
    int ret = solver.ComputeConcentration(x, newton_tol_, newton_maxits_);
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "CALPHADFreeEnergyFunctionsTernary::"
                     "computePhaseConcentrations() "
                     "failed for conc0="
                  << conc0 << ", conc1=" << conc1 << ", hphi=" << hphi
                  << std::endl;
        abort();
    }

    return ret;
}

//-----------------------------------------------------------------------

void CALPHADFreeEnergyFunctionsTernary::energyVsPhiAndC(
    const double temperature, const double* const ceq, const bool found_ceq,
    const double phi_well_scale, const int npts_phi, const int npts_c)
{
    std::clog << "CALPHADFreeEnergyFunctionsTernary::energyVsPhiAndC()..."
              << std::endl;

    const double* const ceqL = &ceq[0];
    const double* const ceqS = &ceq[2];
    // std::clog<<"Input:
    // "<<ceq[0]<<","<<ceq[1]<<","<<ceq[2]<<","<<ceq[3]<<endl;

    double slopec = 0.;
    double fc0    = 0.;
    double fc1    = 0.;
    if (found_ceq)
    {
        // compute slope of f between equilibrium concentrations
        // to add slopec*conc to energy later on

        fc0    = computeFreeEnergy(temperature, &ceqL[0], PhaseIndex::phaseL);
        fc1    = computeFreeEnergy(temperature, &ceqS[0], PhaseIndex::phaseA);
        slopec = -(fc1 - fc0) / (ceqL[1] - ceqL[0]);
    }
    std::clog << std::setprecision(8) << "fc0: " << fc0 << "..."
              << ", fc1: " << fc1 << "..." << std::endl;
    std::clog << "CALPHADFreeEnergyFunctionsTernary: Use slope: " << slopec
              << "..." << std::endl;

    {
        // reset cmin, cmax, deltac
        double c0min = std::min(ceqL[0], ceqS[0]);
        double c0max = std::max(ceqL[0], ceqS[0]);

        double dc0     = c0max - c0min;
        c0min          = std::max(0.25 * c0min, c0min - 0.25 * dc0);
        c0max          = std::min(1. - 0.25 * (1. - c0max), c0max + 0.25 * dc0);
        c0max          = std::max(c0max, c0min + dc0);
        double deltac0 = (c0max - c0min) / (npts_c - 1);

        double c1min = std::min(ceqL[1], ceqS[1]);
        double c1max = std::max(ceqL[1], ceqS[1]);

        double dc1     = c1max - c1min;
        c1min          = std::max(0.25 * c1min, c1min - 0.25 * dc1);
        c1max          = std::min(1. - 0.25 * (1. - c1max), c1max + 0.25 * dc1);
        c1max          = std::max(c1max, c1min + dc1);
        double deltac1 = (c1max - c1min) / (npts_c - 1);

        std::clog << "Range for c0: " << c0min << " to " << c0max << std::endl;
        std::clog << "Range for c1: " << c1min << " to " << c1max << std::endl;

        std::ofstream tfile(fenergy_diag_filename_.data(), std::ios::out);

        printEnergyVsPhiHeader(temperature, npts_phi, npts_c, npts_c, c0min,
            c0max, c1min, c1max, tfile);

        for (int i0 = 0; i0 < npts_c; i0++)
        {
            int i1    = (1. - c0min - deltac0 * i0 - c1min) / deltac1;
            int i1max = i1 < npts_c ? i1 : npts_c;
            for (int i1 = 0; i1 < i1max; i1++)
            {
                double c[2] = { c0min + deltac0 * i0, c1min + deltac1 * i1 };
                printEnergyVsPhi(
                    c, temperature, phi_well_scale, npts_phi, tfile);
            }
        }
    }
}

// Print out free energy as a function of phase
// for given composition and temperature
// File format: ASCII VTK, readble with Visit
void CALPHADFreeEnergyFunctionsTernary::printEnergyVsPhiHeader(
    const double temperature, const int nphi, const int nc0, const int nc1,
    const double c0min, const double c0max, const double c1min,
    const double c1max, std::ostream& os) const
{
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << "Free energy [J/mol] at T=" << temperature << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET STRUCTURED_POINTS" << std::endl;

    os << "DIMENSIONS   " << nphi << " " << nc0 << " " << nc1 << std::endl;
    double asp_ratio_c0 = (nc0 > 1) ? (c0max - c0min) / (nc0 - 1) : 1.;
    double asp_ratio_c1 = (nc1 > 1) ? (c1max - c1min) / (nc1 - 1) : 1.;
    os << "ASPECT_RATIO " << 1. / (nphi - 1) << " " << asp_ratio_c0 << " "
       << asp_ratio_c1 << std::endl;
    os << "ORIGIN        0. " << c0min << " " << c1min << std::endl;
    os << "POINT_DATA   " << nphi * nc0 * nc1 << std::endl;
    os << "SCALARS energy float 1" << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::printEnergyVsPhi(const double* conc,
    const double temperature, const double phi_well_scale, const int npts,
    std::ostream& os)
{
    // std::cout << "CALPHADFreeEnergyFunctionsTernary::printEnergyVsPhi()..."
    // << std::endl;
    const double dphi = 1.0 / (double)(npts - 1);

    // os << "# phi     f(phi)     for c=" << conc
    //           << " and T=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double phi = i * dphi;

        double e       = fchem(phi, conc, temperature);
        const double w = phi_well_scale * well_func(phi);

        os << e + w << std::endl;
    }
    // os << std::endl;
}

//=======================================================================
// compute free energy in [J/mol]
double CALPHADFreeEnergyFunctionsTernary::fchem(
    const double phi, const double* const conc, const double temperature)
{
    const double conc0 = conc[0];
    const double conc1 = conc[1];

    const double hcphi = interp_func(conc_interp_func_type_, phi);

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;
    if ((phi > tol) & (phi < 1. - tol))
    {
        computePhasesFreeEnergies(temperature, hcphi, conc0, conc1, fl, fa);
    }
    else
    {
        // don't solve for phases concentrations, just compute energy
        // in either phase
        double conc[2] = { conc0, conc1 };
        if (phi <= tol)
        {
            fl = computeFreeEnergy(temperature, &conc[0], PhaseIndex::phaseL);
        }
        else
        {
            fa = computeFreeEnergy(temperature, &conc[0], PhaseIndex::phaseA);
        }
    }

    const double hfphi = interp_func(energy_interp_func_type_, phi);

    return (1.0 - hfphi) * fl + hfphi * fa;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsTernary::printEnergyVsComposition(
    const double temperature, const int npts)
{
    const double dc = 1.0 / (double)(npts - 1);

    std::string filename1("Fl");
    filename1 += g_species_phaseL_[0].name();
    filename1 += g_species_phaseL_[2].name();
    filename1 += ".dat";
    std::ofstream os1(filename1.c_str(), std::ios::out);
    os1 << "#phi=0, c1=0, temperature=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        double conc[2];
        conc[0] = i * dc;
        conc[1] = 0.;

        double e = fchem(0., conc, temperature);
        os1 << conc[0] << "\t" << e << std::endl;
    }
    os1 << std::endl;

    std::string filename2("Fs");
    filename2 += g_species_phaseA_[0].name();
    filename2 += g_species_phaseA_[2].name();
    filename2 += ".dat";
    std::ofstream os2(filename2.c_str(), std::ios::out);
    os2 << "#phi=1, c1=0, temperature=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        double conc[2];
        conc[0] = i * dc;
        conc[1] = 0.;

        double e = fchem(1., conc, temperature);
        os2 << conc[0] << "\t" << e << std::endl;
    }
    os2 << std::endl;

    std::string filename3("Fl");
    filename3 += g_species_phaseL_[1].name();
    filename3 += g_species_phaseL_[2].name();
    filename3 += ".dat";
    std::ofstream os3(filename3.c_str(), std::ios::out);
    os3 << "#phi=0, c0=0, temperature=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        double conc[2];
        conc[0] = 0.;
        conc[1] = i * dc;

        double e = fchem(0., conc, temperature);
        os3 << conc[1] << "\t" << e << std::endl;
    }
    os3 << std::endl;

    std::string filename4("Fs");
    filename4 += g_species_phaseA_[1].name();
    filename4 += g_species_phaseA_[2].name();
    filename4 += ".dat";
    std::ofstream os4(filename4.c_str(), std::ios::out);
    os4 << "#phi=1, c0=0, temperature=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        double conc[2];
        conc[0] = 0.;
        conc[1] = i * dc;

        double e = fchem(1., conc, temperature);
        os4 << conc[1] << "\t" << e << std::endl;
    }
    os4 << std::endl;

    std::string filename5("Fl");
    filename5 += g_species_phaseL_[0].name();
    filename5 += g_species_phaseL_[1].name();
    filename5 += ".dat";
    std::ofstream os5(filename5.c_str(), std::ios::out);
    os5 << "#phi=0, temperature=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        double conc[2];
        conc[0] = i * dc;
        conc[1] = 1. - i * dc;

        double e = fchem(0., conc, temperature);
        os5 << conc[0] << "\t" << e << std::endl;
    }
    os5 << std::endl;

    std::string filename6("Fs");
    filename6 += g_species_phaseA_[0].name();
    filename6 += g_species_phaseA_[1].name();
    filename6 += ".dat";
    std::ofstream os6(filename6.c_str(), std::ios::out);
    os6 << "#phi=1, temperature=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        double conc[2];
        conc[0] = i * dc;
        conc[1] = 1. - i * dc;

        double e = fchem(1., conc, temperature);
        os6 << conc[0] << "\t" << e << std::endl;
    }
    os6 << std::endl;
}

void CALPHADFreeEnergyFunctionsTernary::preRunDiagnostics(
    const double T0, const double T1)
{
    std::ofstream os1("FlC0vsT.dat", std::ios::out);
    os1 << "#Species 0, Phase L" << std::endl;
    g_species_phaseL_[0].plotFofT(os1, T0, T1);

    std::ofstream os2("FlC1vsT.dat", std::ios::out);
    os2 << "#Species 1, Phase L" << std::endl;
    g_species_phaseL_[1].plotFofT(os2, T0, T1);

    std::ofstream os3("FlC2vsT.dat", std::ios::out);
    os3 << "#Species 2, Phase L" << std::endl;
    g_species_phaseL_[2].plotFofT(os3, T0, T1);

    std::ofstream os4("FsC0vsT.dat", std::ios::out);
    os4 << "#Species 0, Phase A" << std::endl;
    g_species_phaseA_[0].plotFofT(os4, T0, T1);

    std::ofstream os5("FsC1vsT.dat", std::ios::out);
    os5 << "#Species 1, Phase A" << std::endl;
    g_species_phaseA_[1].plotFofT(os5, T0, T1);

    std::ofstream os6("FsC2vsT.dat", std::ios::out);
    os6 << "#Species 2, Phase A" << std::endl;
    g_species_phaseA_[2].plotFofT(os6, T0, T1);
}
}
