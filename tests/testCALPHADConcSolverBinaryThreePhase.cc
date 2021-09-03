#define CATCH_CONFIG_MAIN

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

TEST_CASE("CALPHAD conc solver binary three phase KKS, two-phase consistancy",
    "[conc solver binary three phase kks, two-phase consistancy]")
{
    // Calculate the inputs and the reference solution
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    double temperature = 1450.;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphadAuNi.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // Get the CALPHAD parameters
    CalphadDataType fA[2];
    CalphadDataType fB[2];
    CalphadDataType Lmix_L[4];
    CalphadDataType Lmix_A[4];

    cafe.computeTdependentParameters(temperature, Lmix_L, Lmix_A, fA, fB);

    // initial guesses
    double c_init0 = 0.5;
    double c_init1 = 0.5;

    double sol_reference[2] = { c_init0, c_init1 };

    // compute concentrations satisfying KKS equations
    const double conc = 0.3;
    const double phi  = 0.5;
    cafe.computePhaseConcentrations(
        temperature, &conc, &phi, &sol_reference[0]);

    // Inputs to the solver
    const double RTinv
        = 1.0 / (Thermo4PFM::gas_constant_R_JpKpmol * temperature);

    double hphi0 = interp_func(conc_interp_func_type, 1.0 - phi);
    double hphi1 = interp_func(conc_interp_func_type, phi);
    double hphi2 = 0.0;

    CalphadDataType Lmix_S0[4];
    CalphadDataType Lmix_S1[4];
    for (int i = 0; i < 4; ++i)
    {
        Lmix_S0[i] = Lmix_A[i];
        Lmix_S1[i] = Lmix_A[i];
    }

    const double tol    = 1.e-8;
    const double alpha  = 0.5; // Using alpha=1 can lead to convergence issues
    const int max_iters = 100;

    // Create and set up the solver
    CalphadDataType fA_threePhase[3];
    CalphadDataType fB_threePhase[3];
    fA_threePhase[0] = fA[0];
    fA_threePhase[1] = fA[1];
    fA_threePhase[2] = fA[1];
    fB_threePhase[0] = fB[0];
    fB_threePhase[1] = fB[1];
    fB_threePhase[2] = fB[1];

    Thermo4PFM::CALPHADConcSolverBinaryThreePhase solver;
    solver.setup(conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_S0, Lmix_S1,
        fA_threePhase, fB_threePhase);

    // Run the solver
    double sol_test[3] = { c_init0, c_init1, 1.0e-8 };
    std::cout << "First test..." << std::endl;
    int ret = solver.ComputeConcentration(sol_test, tol, max_iters, alpha);
    std::cout << "...completed" << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    double residual[3];
    solver.RHS(sol_test, residual);

    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);
    REQUIRE(std::abs(residual[2]) < 1.1 * tol);

    // Make sure that it matches the binary case
    REQUIRE(std::abs(sol_test[0] - sol_reference[0]) < 1.1 * tol);
    REQUIRE(std::abs(sol_test[1] - sol_reference[1]) < 1.1 * tol);

    // ----------
    // Now do the same but for the second solid phase
    hphi0 = interp_func(conc_interp_func_type, 1.0 - phi);
    hphi1 = 0.0;
    hphi2 = interp_func(conc_interp_func_type, phi);

    Thermo4PFM::CALPHADConcSolverBinaryThreePhase solver2;
    solver2.setup(conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_S0, Lmix_S1,
        fA_threePhase, fB_threePhase);

    // Run the solver
    sol_test[0] = c_init0;
    sol_test[1] = 1.0e-8;
    sol_test[2] = c_init1;

    residual[0] = 0.;
    residual[1] = 0.;
    residual[2] = 0.;

    std::cout << "Second test..." << std::endl;
    ret = solver2.ComputeConcentration(sol_test, tol, max_iters, alpha);
    std::cout << "...completed" << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    solver2.RHS(sol_test, residual);

    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);
    REQUIRE(std::abs(residual[2]) < 1.1 * tol);

    // Make sure that it matches the binary case
    REQUIRE(std::abs(sol_test[0] - sol_reference[0]) < 2.0 * tol);
    REQUIRE(std::abs(sol_test[2] - sol_reference[1]) < 2.0 * tol);
}

TEST_CASE("CALPHAD conc solver binary three phase KKS, three-phase convergence",
    "[conc solver binary three phase kks, three-phase convergence]")
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

    Thermo4PFM::CALPHADConcSolverBinaryThreePhase solver;
    solver.setup(
        conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB);

    // Run the solver
    double sol_test[3] = { c_init0, c_init1, c_init2 };
    std::cout << "Third test..." << std::endl;
    int ret = solver.ComputeConcentration(sol_test, tol, max_iters, alpha);
    std::cout << "...completed" << std::endl;
    std::cout << sol_test[0] << " " << sol_test[1] << " " << sol_test[2]
              << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    double residual[3];
    solver.RHS(sol_test, residual);

    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);
    REQUIRE(std::abs(residual[2]) < 1.1 * tol);
}
