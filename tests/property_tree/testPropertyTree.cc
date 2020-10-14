#define CATCH_CONFIG_MAIN

#include "../catch.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace pt = boost::property_tree;

// Simply check that we are able to read correctly a few
// numbers from JSON file
TEST_CASE("Simple JSON read", "[json read]")
{
    std::string input_filename("test.json");

    // Create a root
    pt::ptree troot;

    // Load the json file in this ptree
    pt::read_json(input_filename, troot);

    std::cout << "Reading " << input_filename << std::endl;

    pt::ptree& spA  = troot.get_child("SpeciesA");
    pt::ptree& liqA = spA.get_child("PhaseL");

    // Read values
    std::vector<double> Tcs;
    for (pt::ptree::value_type& tc : liqA.get_child("Tc"))
    {
        Tcs.push_back(tc.second.get_value<double>());
    }

    for (auto tc : Tcs)
        std::cout << " " << tc;
    std::cout << std::endl;

    REQUIRE(Tcs.size() == 5);

    REQUIRE(Tcs[0] == Approx(298.15).epsilon(1.e-12));

    REQUIRE(Tcs[4] == Approx(3000.0).epsilon(1.e-12));
}
