#define CATCH_CONFIG_MAIN

#include "InterpolationType.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("Quadratic binary KKS", "[binary kks]")
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    double Tref = 1693.47;
    double Al   = 143742.0;
    double Aa   = 170630.0;
    double ceql = 0.174;
    double ceqa = 0.154;
    double ml   = -0.0019;
    double ma   = -0.0016;

    Thermo4PFM::QuadraticFreeEnergyFunctionsBinary cafe(Tref, Al, ceql, ml, Aa,
        ceqa, ma, energy_interp_func_type, conc_interp_func_type);

    double sol[2] = { 0., 0. };

    // compute concentrations satisfying KKS equations
    double temperature = 1450.;
    double conc        = 0.3;
    double phi         = 0.4;
    cafe.computePhaseConcentrations(temperature, &conc, &phi, &sol[0]);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cS = " << sol[1] << std::endl;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS);
    std::cout << "   dfS/dcS = " << derivS << std::endl;

    REQUIRE(derivS == Approx(derivL).margin(1.e-5));
}
