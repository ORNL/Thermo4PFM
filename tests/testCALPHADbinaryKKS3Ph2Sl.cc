#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD binary three phase, two sublattice KKS",
    "[binary three phase, two sublattice kks]")
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
            "../thermodynamic_data/calphadAlCuLFccTheta.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    pt::ptree newton_db;
    newton_db.put("alpha", 0.5);
    newton_db.put("max_its", 100);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.8;
    double c_init1 = 0.6;
    double c_init2 = 0.67;

    double sol[3] = { c_init0, c_init1, c_init2 };

    // compute concentrations satisfying KKS equations
    // I think the system at 820K is best behaved in the conc = 0.67-0.8 range
    double conc = 0.7;

    // This phi doesn't quite give Sum(h(phi)) = 1, but it is close
    double phi[3]        = { 0.3, 0.45, 0.5 };
    const double hphi[3] = { interp_func(conc_interp_func_type, phi[0]),
        interp_func(conc_interp_func_type, phi[1]),
        interp_func(conc_interp_func_type, phi[2]) };

    int nit = cafe.computePhaseConcentrations(temperature, &conc, hphi, sol);
    CHECK(nit >= 0);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi[0] << " "
              << phi[1] << " " << phi[2] << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cA = " << sol[1] << std::endl;
    std::cout << "   cB = " << sol[2] << std::endl;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;
    const Thermo4PFM::PhaseIndex pi2 = Thermo4PFM::PhaseIndex::phaseB;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS1;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS1);
    std::cout << "   dfS1/dcS1 = " << derivS1 << std::endl;

    REQUIRE(derivS1 == Approx(derivL).margin(1.e-5));

    double derivS2;
    cafe.computeDerivFreeEnergy(temperature, &sol[2], pi2, &derivS2);
    std::cout << "   dfS2/dcS2 = " << derivS2 << std::endl;

    REQUIRE(derivS2 == Approx(derivL).margin(1.e-5));

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
        REQUIRE(derivFDA == Approx(derivS1).margin(1.e-5));
    }
    {
        double feB      = cafe.computeFreeEnergy(temperature, &sol[2], pi2);
        double ceps     = sol[2] + epsilon;
        double feBeps   = cafe.computeFreeEnergy(temperature, &ceps, pi2);
        double derivFDB = (feBeps - feB) / epsilon;
        std::cout << "FD: dfB/dcB = " << derivFDB << std::endl;
        REQUIRE(derivFDB == Approx(derivS2).margin(1.e-5));
    }

    // check second derivative
    {
        double deriv2L;
        cafe.computeSecondDerivativeFreeEnergy(
            temperature, &sol[0], pi0, &deriv2L);
        std::cout << "d2fL/dcL2 = " << deriv2L << std::endl;
        double derivL0;
        cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL0);
        double derivLp;
        double conc = sol[0] + epsilon;
        cafe.computeDerivFreeEnergy(temperature, &conc, pi0, &derivLp);
        double deriv2FD = (derivLp - derivL0) / epsilon;
        std::cout << "FD: d2fL/dcL2 = " << deriv2FD << std::endl;
        REQUIRE(deriv2FD == Approx(deriv2L).margin(1.e-5));
    }
    {
        double deriv2A;
        cafe.computeSecondDerivativeFreeEnergy(
            temperature, &sol[1], pi1, &deriv2A);
        std::cout << "d2fA/dcA2 = " << deriv2A << std::endl;
        double derivA0;
        cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivA0);
        double derivAp;
        double conc = sol[1] + epsilon;
        cafe.computeDerivFreeEnergy(temperature, &conc, pi1, &derivAp);
        double deriv2FD = (derivAp - derivA0) / epsilon;
        std::cout << "FD: d2fA/dcA2 = " << deriv2FD << std::endl;
        REQUIRE(deriv2FD == Approx(deriv2A).margin(1.e-5));
    }
    {
        double deriv2B;
        cafe.computeSecondDerivativeFreeEnergy(
            temperature, &sol[2], pi2, &deriv2B);
        std::cout << "d2fB/dcB2 = " << deriv2B << std::endl;
        double derivB0;
        cafe.computeDerivFreeEnergy(temperature, &sol[2], pi2, &derivB0);
        double derivBp;
        double conc = sol[2] + epsilon;
        cafe.computeDerivFreeEnergy(temperature, &conc, pi2, &derivBp);
        double deriv2FD = (derivBp - derivB0) / epsilon;
        std::cout << "FD: d2fB/dcB2 = " << deriv2FD << std::endl;
        REQUIRE(deriv2FD == Approx(deriv2B).margin(1.e-5));
    }
}

TEST_CASE("CALPHAD binary three phase, two sublattice KKS #2",
    "[binary three phase, two sublattice kks #2]")
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
            "../thermodynamic_data/calphadAlCuLFccTheta.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    pt::ptree newton_db;
    newton_db.put("alpha", 1.0);
    newton_db.put("max_its", 20000);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.79267;
    double c_init1 = 0.79267;
    double c_init2 = 0.79267;

    double sol[3] = { c_init0, c_init1, c_init2 };

    // compute concentrations satisfying KKS equations
    // I think the system at 820K is best behaved in the conc = 0.67-0.8 range
    double conc = 0.79267;

    // This phi doesn't quite give Sum(h(phi)) = 1, but it is close
    double phi[3] = { 0.938956, 5.88418e-15, 0.0610444 };

    cafe.computePhaseConcentrations(temperature, &conc, phi, sol);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi[0] << " "
              << phi[1] << " " << phi[2] << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cA = " << sol[1] << std::endl;
    std::cout << "   cB = " << sol[2] << std::endl;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;
    const Thermo4PFM::PhaseIndex pi2 = Thermo4PFM::PhaseIndex::phaseB;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS1;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS1);
    std::cout << "   dfS1/dcS1 = " << derivS1 << std::endl;

    // REQUIRE(derivS1 == Approx(derivL).margin(1.e-5));

    double derivS2;
    cafe.computeDerivFreeEnergy(temperature, &sol[2], pi2, &derivS2);
    std::cout << "   dfS2/dcS2 = " << derivS2 << std::endl;

    REQUIRE(derivS2 == Approx(derivL).margin(1.e-5));
}

TEST_CASE("CALPHAD binary three phase, two sublattice KKS #3",
    "[binary three phase, two sublattice kks #3]")
{
    // This case may require a restart of the solver with new initial conditions
    // to pass.

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

    pt::ptree newton_db;
    newton_db.put("alpha", 1.0);
    newton_db.put("max_its", 30);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.800000011921;
    double c_init1 = 0.800000011921;
    double c_init2 = 0.800000011921;

    double sol[3] = { c_init0, c_init1, c_init2 };

    // compute concentrations satisfying KKS equations
    // I think the system at 820K is best behaved in the conc = 0.67-0.8 range
    double conc = 0.800000011921;

    // Want phi that corresponds to hphi = [0.982551, 0.0174491, 4.79499e-25]
    // Want phi that corresponds to hphi = [1.0, 6.18271e-09, 3.69779e-32]
    // double phi[3] = { 0.998, 8.5e-4, 1.5e-11 };
    double phi[3] = { 0.998, 8.5e-4, 1.5e-11 };

    cafe.computePhaseConcentrations(temperature, &conc, phi, sol);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi[0] << " "
              << phi[1] << " " << phi[2] << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cA = " << sol[1] << std::endl;
    std::cout << "   cB = " << sol[2] << std::endl;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;
    const Thermo4PFM::PhaseIndex pi2 = Thermo4PFM::PhaseIndex::phaseB;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS1;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS1);
    std::cout << "   dfS1/dcS1 = " << derivS1 << std::endl;

    REQUIRE(derivS1 == Approx(derivL).margin(1.e-5));

    double derivS2;
    cafe.computeDerivFreeEnergy(temperature, &sol[2], pi2, &derivS2);
    std::cout << "   dfS2/dcS2 = " << derivS2 << std::endl;

    REQUIRE(derivS2 == Approx(derivL).margin(1.e-5));
}
