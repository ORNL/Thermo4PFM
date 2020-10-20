#define CATCH_CONFIG_MAIN

#include "CALPHADSpeciesPhaseGibbsEnergy.h"

#include "catch.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace pt = boost::property_tree;

// Read CALPHAD database for Ni
// Compute a few energies as a function of temperature
// near melting temperature
TEST_CASE("CALPHAD Gibbs energy", "[gibbs energy]")
{
    std::cout << std::setprecision(12);

    // Ni melting temperature
    double temperature = 1728.;

    // read calphad database
    std::string calphad_filename("../thermodynamic_data/calphadAuNi.json");
    pt::ptree calphad_db;
    std::cout << "Read " << calphad_filename << std::endl;
    try
    {
        pt::read_json(calphad_filename, calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    std::cout << "File content:" << std::endl;
    for (auto& v : calphad_db)
        std::cout << v.first << std::endl;

    // Ni is species "B" in calphadAuNi.json
    std::cout << "Read species B data" << std::endl;
    pt::ptree spB;
    try
    {
        spB = calphad_db.get_child("SpeciesB");
    }
    catch (std::exception& e)
    {
        std::cerr << "Read species: exception caught: " << e.what()
                  << std::endl;
    }

    pt::ptree liqB;
    try
    {
        liqB = spB.get_child("PhaseL");
    }
    catch (std::exception& e)
    {
        std::cerr << "Read phase: exception caught: " << e.what() << std::endl;
    }

    pt::ptree solB;
    try
    {
        solB = spB.get_child("PhaseA");
    }
    catch (std::exception& e)
    {
        std::cerr << "Read phase: exception caught: " << e.what() << std::endl;
    }

    // liquid
    CALPHADSpeciesPhaseGibbsEnergy energyNiLiquid;
    energyNiLiquid.initialize("NiL", liqB);

    double el = energyNiLiquid.fenergy(temperature);
    std::cout << "Ni liquid, Energy at T = " << temperature << " : " << el
              << std::endl;

    // solid
    CALPHADSpeciesPhaseGibbsEnergy energyNiSolid;
    energyNiSolid.initialize("NiS", solB);

    double es = energyNiSolid.fenergy(temperature);
    std::cout << "Ni solid, Energy at T = " << temperature << " : " << es
              << std::endl;

    // at the meting temperature, energies in liquid and solid should be equal
    CHECK(es == Approx(el).margin(0.2));

    double belowT = temperature - 1.;
    el            = energyNiLiquid.fenergy(belowT);
    std::cout << "Ni liquid, Energy at T = " << belowT << " : " << el
              << std::endl;
    es = energyNiSolid.fenergy(belowT);
    std::cout << "Ni solid, Energy at T = " << belowT << " : " << es
              << std::endl;

    // below  the meting temperature, energies in solid should be lower than
    // in liquid
    CHECK(es < el - 10.);

    double aboveT = temperature + 1.;
    el            = energyNiLiquid.fenergy(aboveT);
    std::cout << "Ni liquid, Energy at T = " << aboveT << " : " << el
              << std::endl;
    es = energyNiSolid.fenergy(aboveT);
    std::cout << "Ni solid, Energy at T = " << aboveT << " : " << es
              << std::endl;
    CHECK(el < es - 8.);
}
