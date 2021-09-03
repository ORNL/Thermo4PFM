#define CATCH_CONFIG_MAIN

#include "CALPHADConcSolverBinary3Ph2Sl.h"
#include "CALPHADConcSolverBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"

#include "InterpolationType.h"
#include "PhysicalConstants.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD conc solver binary 3 phase, 2 sublattice KKS, "
          "single-sublattice consistency",
    "[conc solver binary 3 phase, 2 sublattice kks, single-sublattice "
    "consistency]")
{
    // Calculate the inputs and the reference solution
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    double temperature = 900.;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json(
            "../thermodynamic_data/calphadAlCuLFccBcc.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    // First calculate the reference solution

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // Get the CALPHAD parameters
    CalphadDataType fA[3];
    CalphadDataType fB[3];
    CalphadDataType Lmix_L[4];
    CalphadDataType Lmix_A[4];
    CalphadDataType Lmix_B[4];

    cafe.computeTdependentParameters(
        temperature, Lmix_L, Lmix_A, Lmix_B, fA, fB);

    // initial guesses
    double c_init0 = 0.9;
    double c_init1 = 0.9;
    double c_init2 = 0.9;

    // compute concentrations satisfying KKS equations
    const double conc = 0.9;

    // Inputs to the solver
    const double RTinv
        = 1.0 / (Thermo4PFM::gas_constant_R_JpKpmol * temperature);

    double hphi0 = interp_func(conc_interp_func_type, 0.5);
    double hphi1 = interp_func(conc_interp_func_type, 0.4);
    double hphi2 = interp_func(conc_interp_func_type, 0.1);

    const double tol    = 1.e-8;
    const double alpha  = 0.1; // Using alpha=1 can lead to convergence issues
    const int max_iters = 10000;

    // Create and set up the solver

    Thermo4PFM::CALPHADConcSolverBinaryThreePhase solver_ref;
    solver_ref.setup(
        conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB);

    // Run the solver
    double sol_ref[3] = { c_init0, c_init1, c_init2 };
    int ret = solver_ref.ComputeConcentration(sol_ref, tol, max_iters, alpha);
    REQUIRE(ret >= 0);

    // Now do the same calculation with the two-sublattice solver
    int p[3];
    int q[3];
    for (int i = 0; i < 3; ++i)
    {
        p[i] = 0;
        q[i] = 1;
    }
    Thermo4PFM::CALPHADConcSolverBinary3Ph2Sl solver;
    solver.setup(
        conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB, p, q);

    // Run the solver
    double sol_test[3] = { c_init0, c_init1, c_init2 };
    ret = solver.ComputeConcentration(sol_test, tol, max_iters, alpha);
    REQUIRE(ret >= 0);

    // Check consistency
    REQUIRE(sol_test[0] == Approx(sol_ref[0]).margin(1.e-5));
    REQUIRE(sol_test[1] == Approx(sol_ref[1]).margin(1.e-5));
    REQUIRE(sol_test[2] == Approx(sol_ref[2]).margin(1.e-5));
}
