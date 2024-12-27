#define CATCH_CONFIG_MAIN

#include "ParabolicEqConcSolverBinary.h"

#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <string>

TEST_CASE("Parabolic conc solver binary three phase KKS, two-phase consistancy",
    "[conc solver binary three phase kks, two-phase consistancy]")
{
    // Parameters from Yang, Wang, Xing, Huang, Mat. Today Comm. 28 (2021)
    // 102712 Ti-Al
    double coeffA[3][2] = { { 22130.9, -4.4354 }, { -10192.4, 2.6104 },
        { -10276.25, -8.7833 } };
    double coeffB[3][2]
        = { { 22090.9, -4.4330 }, { -9841.2, 2.4252 }, { -10439.2, -8.8977 } };
    double coeffL[3][2]
        = { { 17183.8, -3.9748 }, { -8659.95, 2.4447 }, { -10439.2, -9.6344 } };
    const double Tref = 1500.;

    const double temperature = 1520.;
    std::cout << "temperature = " << temperature << std::endl;

    // Create and set up the solver for phases L and A
    Thermo4PFM::ParabolicEqConcSolverBinary solver;
    solver.setup(temperature - Tref, coeffL, coeffA);
    const double toln   = 1.e-8;
    const int max_iters = 10;

    // Run the solver
    double sol_test[2] = { 0.5, 0.5 };
    std::cout << "Test L,A..." << std::endl;
    int ret = solver.ComputeConcentration(sol_test, toln, max_iters);
    std::cout << ret << " iterations" << std::endl;
    std::cout << "Solution: " << sol_test[0] << " " << sol_test[1] << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    double residual[2];
    solver.RHS(sol_test, residual);

    const double tol = 1.e-8;
    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);

    // ----------
    // Now do the same but for the second solid phase
    solver.setup(temperature - Tref, coeffL, coeffB);

    // Run the solver
    sol_test[0] = 0.5;
    sol_test[1] = 0.5;
    std::cout << "Test L,B..." << std::endl;
    ret = solver.ComputeConcentration(sol_test, toln, max_iters);
    std::cout << ret << " iterations" << std::endl;
    std::cout << "Solution: " << sol_test[0] << " " << sol_test[1] << std::endl;
    REQUIRE(ret >= 0);

    residual[0] = 0.;
    residual[1] = 0.;
    solver.RHS(sol_test, residual);

    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);

    // Inputs to the solver
    // setup the solver
    solver.setup(temperature - Tref, coeffA, coeffB);

    // Run the solver
    sol_test[0] = 0.5;
    sol_test[1] = 0.5;
    std::cout << "Test A,B..." << std::endl;
    ret = solver.ComputeConcentration(sol_test, toln, max_iters);
    std::cout << ret << " iterations" << std::endl;
    std::cout << "Solution: " << sol_test[0] << " " << sol_test[1] << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the RHS
    residual[0] = 0.;
    residual[1] = 0.;
    solver.RHS(sol_test, residual);

    REQUIRE(std::abs(residual[0]) < 1.1 * tol);
    REQUIRE(std::abs(residual[1]) < 1.1 * tol);
}
