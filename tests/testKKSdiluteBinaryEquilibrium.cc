#define CATCH_CONFIG_MAIN

#include "KKSFreeEnergyFunctionDiluteBinary.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

namespace pt = boost::property_tree;

TEST_CASE("Dilute binary equilibrium", "[dilute binary equilibrium]")
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    // read dilute alloy informaton
    pt::ptree conc_db;
    try
    {
        pt::read_json("dilute_binary.json", conc_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }
    pt::write_json(std::clog, conc_db);

    double temperature = 1710.;

    Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary cafe(
        conc_db, energy_interp_func_type, conc_interp_func_type);

    const double tol = 1.e-5;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;

    // compute equilibrium compositions
    double ceq[2];
    cafe.computeCeqT(temperature, pi0, pi1, &ceq[0]);
    std::cout << "   ceL = " << ceq[0] << std::endl;
    std::cout << "   ceS = " << ceq[1] << std::endl;

    // test if chemical potentials are equal for equilibrium compositions
    double derivL;
    double derivS;
    cafe.computeDerivFreeEnergy(temperature, &ceq[0], pi0, &derivL);
    cafe.computeDerivFreeEnergy(temperature, &ceq[1], pi1, &derivS);
    std::cout << "   dfL/dcL = " << derivL << std::endl;
    std::cout << "   dfS/dcS = " << derivS << std::endl;
    CHECK(derivS == Approx(derivL).margin(tol));

    // test if driving force is 0 for equilibrium compositions
    double fl = cafe.computeFreeEnergy(temperature, &ceq[0], pi0, false);
    double fa = cafe.computeFreeEnergy(temperature, &ceq[1], pi1, false);
    std::cout << "   fL = " << fl << std::endl;
    std::cout << "   fS = " << fa << std::endl;
    double diff = fa - fl - derivS * (ceq[1] - ceq[0]);
    CHECK(std::abs(diff) == Approx(0.).margin(tol));

    // compute second derivatives for info only
    std::vector<double> d2fdc2(1);
    cafe.computeSecondDerivativeFreeEnergy(
        temperature, &ceq[0], pi0, d2fdc2.data());
    std::cout << "-------------------------------" << std::endl;
    std::cout << "Second derivatives" << std::endl;
    std::cout << "At ceL: " << d2fdc2[0] << std::endl;
    cafe.computeSecondDerivativeFreeEnergy(
        temperature, &ceq[1], pi1, d2fdc2.data());
    std::cout << "At ceS: " << d2fdc2[0] << std::endl;
}
