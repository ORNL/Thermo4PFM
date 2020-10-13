#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "USAGE:  " << argv[0] << " <input filename> ";
        return (-1);
    }

    std::string input_filename = argv[1];

    // Create a root
    pt::ptree troot;

    // Load the json file in this ptree
    pt::read_json(input_filename, troot);

    std::cout << "Reading test.json..." << std::endl;

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
}
