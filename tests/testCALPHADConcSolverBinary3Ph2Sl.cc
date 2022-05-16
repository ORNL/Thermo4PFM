#define CATCH_CONFIG_MAIN

#include "CALPHADConcSolverBinary3Ph2Sl.h"
#include "CALPHADConcSolverBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
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
    const double RT    = GASCONSTANT_R_JPKPMOL * temperature;
    const double RTinv = 1.0 / RT;

    double hphi0 = interp_func(conc_interp_func_type, 0.5);
    double hphi1 = interp_func(conc_interp_func_type, 0.4);
    double hphi2 = interp_func(conc_interp_func_type, 0.1);

    const double tol    = 1.e-8;
    const double alpha  = 0.1; // Using alpha=1 can lead to convergence issues
    const int max_iters = 10000;

    // Create and set up the solver

    Thermo4PFM::CALPHADConcSolverBinaryThreePhase solver_ref;
    solver_ref.setup(
        conc, hphi0, hphi1, hphi2, RT, Lmix_L, Lmix_A, Lmix_B, fA, fB);

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

TEST_CASE("CALPHAD conc solver binary 3 phase, 2 sublattice KKS, "
          "convergence",
    "[conc solver binary 3 phase, 2 sublattice kks, convergence]")
{
    // Calculate the inputs and the reference solution
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    double temperature = 820.;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json(
            "../thermodynamic_data/calphadAlCuLFccTheta.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    // Test 1

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
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
    double c_init0 = 0.8;
    double c_init1 = 0.6;
    double c_init2 = 0.67;

    // compute concentrations satisfying KKS equations

    // I think the system at 820K is best behaved in the conc = 0.67-0.8 range
    double conc = 0.7;

    // Inputs to the solver
    const double RTinv = 1.0 / (GASCONSTANT_R_JPKPMOL * temperature);

    // NOTE: The sum of hphi should equal one, which is not necessarily true for
    // phi
    double hphi0 = 0.1;
    double hphi1 = 0.4;
    double hphi2 = 0.5;

    double tol    = 1.e-8;
    double alpha  = 0.5; // Using alpha=1 can lead to convergence issues
    int max_iters = 100;

    int p[3];
    int q[3];
    for (int i = 0; i < 3; ++i)
    {
        p[i] = 0;
        q[i] = 1;
    }

    p[2] = 2;

    Thermo4PFM::CALPHADConcSolverBinary3Ph2Sl solver;
    solver.setup(
        conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB, p, q);

    // Run the solver
    double sol_test[3] = { c_init0, c_init1, c_init2 };
    int ret = solver.ComputeConcentration(sol_test, tol, max_iters, alpha);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << std::endl;
    std::cout << "   cL = " << sol_test[0] << std::endl;
    std::cout << "   cA = " << sol_test[1] << std::endl;
    std::cout << "   cB = " << sol_test[2] << std::endl;

    REQUIRE(ret >= 0);

    // Test 2
    max_iters = 20000;
    alpha     = 1.0;

    sol_test[0] = 0.79267;
    sol_test[1] = 0.79267;
    sol_test[2] = 0.79267;

    hphi0 = 0.989276;
    hphi1 = 1.03871e-28;
    hphi2 = 0.0107243;

    conc = 0.79267;

    solver.setup(
        conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB, p, q);

    // Run the solver
    ret = solver.ComputeConcentration(sol_test, tol, max_iters, alpha);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << std::endl;
    std::cout << "   cL = " << sol_test[0] << std::endl;
    std::cout << "   cA = " << sol_test[1] << std::endl;
    std::cout << "   cB = " << sol_test[2] << std::endl;
    std::cout << "Newton iterations = " << ret << std::endl;
    REQUIRE(ret >= 0);

    // Test 3

    max_iters = 2000;
    alpha     = 1.0;

    sol_test[0] = 0.7;
    sol_test[1] = 0.9;
    sol_test[2] = 0.794081;

    hphi0 = 0.982551;
    hphi1 = 0.0174491;
    hphi2 = 4.79499e-25;

    conc = 0.794081;

    solver.setup(
        conc, hphi0, hphi1, hphi2, RTinv, Lmix_L, Lmix_A, Lmix_B, fA, fB, p, q);

    // Run the solver
    ret = solver.ComputeConcentration(sol_test, tol, max_iters, alpha);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << std::endl;
    std::cout << "   cL = " << sol_test[0] << std::endl;
    std::cout << "   cA = " << sol_test[1] << std::endl;
    std::cout << "   cB = " << sol_test[2] << std::endl;
    std::cout << "Newton iterations = " << ret << std::endl;
    REQUIRE(ret >= 0);
}
