#define CATCH_CONFIG_MAIN

#include "Phases.h"
#include "QuadraticFreeEnergyFunctionsBinaryThreePhase.h"

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

    // Inputs to the solver
    double hphi[3] = { 0.2, 0.3, 0.5 };

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    // Create and set up the solver
    Thermo4PFM::QuadraticFreeEnergyFunctionsBinaryThreePhase qfe(Al, ceql, Aa,
        ceqa, Ab, ceqb, energy_interp_func_type, conc_interp_func_type);

    // Run the solver
    double sol_test[3] = { 0., 0., 0. };
    std::cout << "First test..." << std::endl;
    int ret = qfe.computePhaseConcentrations(0., &conc, hphi, sol_test);
    std::cout << sol_test[0] << " " << sol_test[1] << " " << sol_test[2]
              << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the derivative computation
    Thermo4PFM::PhaseIndex pl = Thermo4PFM::PhaseIndex::phaseL;
    double derivl;
    qfe.computeDerivFreeEnergy(0., sol_test, pl, &derivl);
    std::cout << "dfl/dcl = " << derivl << std::endl;

    Thermo4PFM::PhaseIndex pa = Thermo4PFM::PhaseIndex::phaseA;
    double deriva;
    qfe.computeDerivFreeEnergy(0., sol_test + 1, pa, &deriva);
    std::cout << "dfa/dca = " << deriva << std::endl;

    Thermo4PFM::PhaseIndex pb = Thermo4PFM::PhaseIndex::phaseB;
    double derivb;
    qfe.computeDerivFreeEnergy(0., sol_test + 2, pb, &derivb);
    std::cout << "dfb/dcb = " << derivb << std::endl;

    const double tol = 1.e-8;
    REQUIRE(std::abs(deriva - derivl) < tol);
    REQUIRE(std::abs(derivb - derivl) < tol);
}
