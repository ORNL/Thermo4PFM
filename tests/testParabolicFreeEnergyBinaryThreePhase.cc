#define CATCH_CONFIG_MAIN

#include "ParabolicFreeEnergyFunctionsBinaryThreePhase.h"
#include "Phases.h"

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
        = { { 22090.9, -4.4330 }, { -9841.2, 2.4252 }, { -10430.25, -8.8977 } };
    double coeffL[3][2]
        = { { 17183.8, -3.9748 }, { -8659.95, 2.4447 }, { -10439.2, -9.6344 } };
    const double Tref = 1500.;

    const double conc = 0.46;
    std::cout << "conc = " << conc << std::endl;

    const double temperature = 1520.;
    std::cout << "temperature = " << temperature << std::endl;

    // Inputs to the solver
    double hphi[3] = { 0.32, 0.3, 0.38 };

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    // Create and set up the solver
    Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase qfe(Tref, coeffL,
        coeffA, coeffB, energy_interp_func_type, conc_interp_func_type);

    const int npts = 100;
    std::ofstream os("FvsC_parabolic.dat", std::ios::out);
    std::cout << " Compute energies..." << std::endl;
    const double cmin = 0.37;
    const double cmax = 0.55;
    qfe.printEnergyVsComposition(temperature, os, cmin, cmax, npts);

    // Run the solver
    double sol_test[3] = { 0., 0., 0. };
    std::cout << "First test..." << std::endl;
    int nit
        = qfe.computePhaseConcentrations(temperature, &conc, hphi, sol_test);
    std::cout << "nit = " << nit << std::endl;
    std::cout << sol_test[0] << " " << sol_test[1] << " " << sol_test[2]
              << std::endl;
    REQUIRE(nit >= 0);

    // Plug the solution back into the derivative computation
    Thermo4PFM::PhaseIndex pl = Thermo4PFM::PhaseIndex::phaseL;
    double derivl;
    qfe.computeDerivFreeEnergy(temperature, sol_test, pl, &derivl);
    std::cout << "dfl/dcl = " << derivl << std::endl;

    Thermo4PFM::PhaseIndex pa = Thermo4PFM::PhaseIndex::phaseA;
    double deriva;
    qfe.computeDerivFreeEnergy(temperature, sol_test + 1, pa, &deriva);
    std::cout << "dfa/dca = " << deriva << std::endl;

    Thermo4PFM::PhaseIndex pb = Thermo4PFM::PhaseIndex::phaseB;
    double derivb;
    qfe.computeDerivFreeEnergy(temperature, sol_test + 2, pb, &derivb);
    std::cout << "dfb/dcb = " << derivb << std::endl;

    const double tol = 1.e-8;
    REQUIRE(std::abs(deriva - derivl) < tol);
    REQUIRE(std::abs(derivb - derivl) < tol);
}
