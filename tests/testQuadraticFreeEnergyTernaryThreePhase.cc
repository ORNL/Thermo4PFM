#define CATCH_CONFIG_MAIN

#include "Phases.h"
#include "QuadraticFreeEnergyFunctionsTernaryThreePhase.h"

#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <string>

TEST_CASE("Quadratic conc solver binary three phase KKS, two-phase consistancy",
    "[conc solver binary three phase kks, two-phase consistancy]")
{
    const double Al[2] = { 100., 130. };
    const double Aa[2] = { 110., 120. };
    const double Ab[2] = { 80., 90. };

    const double ceql[2] = { 0.1, 0.3 };
    const double ceqa[2] = { 0.3, 0.6 };
    const double ceqb[2] = { 0.8, 0.7 };

    // Inputs to the solver
    const double conc[2] = { 0.3, 0.2 };
    std::cout << "conc = " << conc[0] << ", " << conc[1] << std::endl;

    double hphi[3] = { 0.2, 0.3, 0.5 };

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    // Create and set up the solver
    Thermo4PFM::QuadraticFreeEnergyFunctionsTernaryThreePhase qfe(Al, ceql, Aa,
        ceqa, Ab, ceqb, energy_interp_func_type, conc_interp_func_type);

    // Run the solver
    double sol_test[6];
    int ret = qfe.computePhaseConcentrations(0., conc, hphi, sol_test);
    std::cout << "Nb. iterations: " << ret << std::endl;
    std::cout << "Phase concentrations:" << std::endl;
    std::cout << sol_test[0] << ", " << sol_test[1] << ", " << sol_test[2]
              << ", " << sol_test[3] << ", " << sol_test[4] << ", "
              << sol_test[5] << std::endl;
    REQUIRE(ret >= 0);

    // Plug the solution back into the derivative computation
    Thermo4PFM::PhaseIndex pl = Thermo4PFM::PhaseIndex::phaseL;
    double derivl[2];
    qfe.computeDerivFreeEnergy(0., sol_test, pl, derivl);
    std::cout << "dfl/dcl = " << derivl[0] << ", " << derivl[1] << std::endl;

    Thermo4PFM::PhaseIndex pa = Thermo4PFM::PhaseIndex::phaseA;
    double deriva[2];
    qfe.computeDerivFreeEnergy(0., sol_test + 2, pa, deriva);
    std::cout << "dfa/dca = " << deriva[0] << ", " << deriva[1] << std::endl;

    Thermo4PFM::PhaseIndex pb = Thermo4PFM::PhaseIndex::phaseB;
    double derivb[2];
    qfe.computeDerivFreeEnergy(0., sol_test + 4, pb, derivb);
    std::cout << "dfb/dcb = " << derivb[0] << ", " << derivb[1] << std::endl;

    const double tol = 1.e-5;
    REQUIRE(std::abs(deriva[0] - derivl[0]) < tol);
    REQUIRE(std::abs(derivb[0] - derivl[0]) < tol);
}
