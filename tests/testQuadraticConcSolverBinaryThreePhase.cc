#define CATCH_CONFIG_MAIN

#include "QuadraticConcSolverBinaryThreePhase.h"

#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <string>

TEST_CASE("Quadratic conc solver binary three phase KKS, two-phase consistancy",
    "[conc solver binary three phase kks, two-phase consistancy]")
{
    const double Al = 10000.;
    const double Aa = 10000.;
    const double Ab = 10000.;

    const double ceql = 0.1;
    const double ceqa = 0.3;
    const double ceqb = 0.8;

    const double conc = 0.3;
    std::cout << "conc = " << conc << std::endl;
    const double phi = 0.5;

    // Inputs to the solver
    double hphi0 = 1.0 - phi;
    double hphi1 = phi;
    double hphi2 = 0.0;

    // Create and set up the solver
    Thermo4PFM::QuadraticConcSolverBinaryThreePhase solver;
    solver.setup(conc, hphi0, hphi1, hphi2, Al, ceql, Aa, ceqa, Ab, ceqb);

    // Run the solver
    double sol_test[3] = { 0., 0., 0. };
    std::cout << "First test..." << std::endl;
    int ret = solver.ComputeConcentration(sol_test);
    std::cout << sol_test[0] << " " << sol_test[1] << " " << sol_test[2]
              << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    double residual[3];
    solver.RHS(sol_test, residual);

    const double tol = 1.e-8;
    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);
    REQUIRE(std::abs(residual[2]) < 1.1 * tol);

    // ----------
    // Now do the same but for the second solid phase
    hphi0 = 1.0 - phi;
    hphi1 = 0.0;
    hphi2 = phi;

    Thermo4PFM::QuadraticConcSolverBinaryThreePhase solver2;
    solver2.setup(conc, hphi0, hphi1, hphi2, Al, ceql, Aa, ceqa, Ab, ceqb);

    // Run the solver
    sol_test[0] = 0.;
    sol_test[1] = 0.;
    sol_test[2] = 0.;

    residual[0] = 0.;
    residual[1] = 0.;
    residual[2] = 0.;

    std::cout << "Second test..." << std::endl;
    ret = solver2.ComputeConcentration(sol_test);
    std::cout << sol_test[0] << " " << sol_test[1] << " " << sol_test[2]
              << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    solver2.RHS(sol_test, residual);

    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);
    REQUIRE(std::abs(residual[2]) < 1.1 * tol);

    // Inputs to the solver
    hphi0 = 0.5;
    hphi1 = 0.4;
    hphi2 = 0.1;

    // setup the solver
    solver.setup(conc, hphi0, hphi1, hphi2, Al, ceql, Aa, ceqa, Ab, ceqb);

    // Run the solver
    sol_test[0] = 0.;
    sol_test[1] = 0.;
    sol_test[2] = 0.;
    std::cout << "Third test..." << std::endl;
    ret = solver.ComputeConcentration(sol_test);
    std::cout << sol_test[0] << " " << sol_test[1] << " " << sol_test[2]
              << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    solver.RHS(sol_test, residual);

    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);
    REQUIRE(std::abs(residual[2]) < 1.1 * tol);
}
