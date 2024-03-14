#define CATCH_CONFIG_MAIN

#include "CALPHADConcSolverBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

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

TEST_CASE("CALPHAD conc solver binary three phase KKS stochio",
    "[conc solver binary three phase kks stochio]")
{
    // Calculate the inputs and the reference solution
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::LINEAR;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    // very small temperature to have a very small contribution of
    // pure mixing term
    // Cannot be 0 as it would trigger nans in f(T) evaluations
    double temperature = 1.e-8;

    std::cout << "Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphad3phases.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // Inputs to the solver: 3 phase fractions, 1 composition (binary)
    double hphi[3] = { interp_func(conc_interp_func_type, 0.5),
        interp_func(conc_interp_func_type, 0.4),
        interp_func(conc_interp_func_type, 0.1) };
    double conc    = 0.5;

    const double tol = 1.e-8;

    // set initial guess for iterative solver
    double x[3] = { 0.5, 0.5, 0.5 };

    // compute concentrations satisfying KKS equations
    int nits = cafe.computePhaseConcentrations(temperature, &conc, hphi, x);
    std::cout << "Solution : " << x[0] << " " << x[1] << " " << x[2]
              << std::endl;
    std::cout << "Number of iterations : " << nits << std::endl;
    REQUIRE(nits > 0);

    // check residual of weighted composition equation
    double residual = conc - hphi[0] * x[0] - hphi[1] * x[1] - hphi[2] * x[2];
    std::cout << "Residual : " << residual << std::endl;
    REQUIRE(std::abs(residual) < 1.1 * tol);

    // save solution
    double sol[3] = { x[0], x[1], x[2] };

    //////////////////////////////////////////////////////////////////////////
    std::cout << "Now stochio test..." << std::endl;
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB cafeStochioB(
        x[2], calphad_db, newton_db, energy_interp_func_type,
        conc_interp_func_type);

    // reset initial guess
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = -1.;
    nits = cafeStochioB.computePhaseConcentrations(temperature, &conc, hphi, x);
    std::cout << "StochioB Solution : " << x[0] << " " << x[1] << " " << x[2]
              << std::endl;
    std::cout << "Number of iterations : " << nits << std::endl;

    // check if solver converged
    REQUIRE(nits > 0);

    // check if we find the same solution
    REQUIRE(std::abs(x[0] - sol[0]) < tol);
    REQUIRE(std::abs(x[1] - sol[1]) < tol);
}
