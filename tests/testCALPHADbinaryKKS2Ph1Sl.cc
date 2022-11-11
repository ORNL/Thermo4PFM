#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD binary KKS", "[binary kks]")
{
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
            "../thermodynamic_data/calphadAlCuLTheta.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.5;
    double c_init1 = 0.5;

    double sol[2] = { c_init0, c_init1 };

    // compute concentrations satisfying KKS equations
    double conc = 0.7;
    double phi  = 0.5;
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

    double derivA;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivA);
    std::cout << "   dfA/dcA = " << derivA << std::endl;

    REQUIRE(derivA == Approx(derivL).margin(1.e-5));

    // check derivatives are consistent with energy
    const double epsilon = 1.e-8;
    {
        double feL      = cafe.computeFreeEnergy(temperature, &sol[0], pi0);
        double ceps     = sol[0] + epsilon;
        double feLeps   = cafe.computeFreeEnergy(temperature, &ceps, pi0);
        double derivFDL = (feLeps - feL) / epsilon;
        std::cout << "FD: dfL/dcL = " << derivFDL << std::endl;
        REQUIRE(derivFDL == Approx(derivL).margin(1.e-5));
    }
    {
        double feA      = cafe.computeFreeEnergy(temperature, &sol[1], pi1);
        double ceps     = sol[1] + epsilon;
        double feAeps   = cafe.computeFreeEnergy(temperature, &ceps, pi1);
        double derivFDA = (feAeps - feA) / epsilon;
        std::cout << "FD: dfA/dcA = " << derivFDA << std::endl;
        REQUIRE(derivFDA == Approx(derivA).margin(1.e-5));
    }

    // check 2nd derivatives are consistent with 1st derivatives
    {
        double dL;
        cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &dL);
        double ceps = sol[0] + epsilon;
        double dLeps;
        cafe.computeDerivFreeEnergy(temperature, &ceps, pi0, &dLeps);
        double derivFDL = (dLeps - dL) / epsilon;
        std::cout << "FD: d2fL/dcL2 = " << derivFDL << std::endl;

        double dL2;
        cafe.computeSecondDerivativeFreeEnergy(temperature, &ceps, pi0, &dL2);
        std::cout << "d2fL/dcL2 = " << dL2 << std::endl;
        REQUIRE(derivFDL == Approx(dL2).margin(1.e-5));
    }
    {
        double dA;
        cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &dA);
        double ceps = sol[1] + epsilon;
        double dAeps;
        cafe.computeDerivFreeEnergy(temperature, &ceps, pi1, &dAeps);
        double derivFDA = (dAeps - dA) / epsilon;
        std::cout << "FD: d2fA/dcA2 = " << derivFDA << std::endl;

        double dA2;
        cafe.computeSecondDerivativeFreeEnergy(temperature, &ceps, pi1, &dA2);
        std::cout << "d2fA/dcA2 = " << dA2 << std::endl;
        REQUIRE(derivFDA == Approx(dA2).margin(1.e-5));
    }
}
