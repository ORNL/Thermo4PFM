#include "CALPHADFreeEnergyFunctionsTernary.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    double temperature = 1570;
    std::cout << "Temperature: " << temperature << std::endl;
    double nominalc[2] = { 0.6, 0.035 };

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("calphadFeNbNi_Mathon_et_al.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    const int nphi = 20;

    const double dphi = 1. / (double)nphi;

    // initial guesses
    double init_guess[5] = { 0.5, 0.05, 0.6, 0.01, 0.5 };
    double lceq[5]       = { init_guess[0], init_guess[1], // liquid
        init_guess[2], init_guess[3], // solid
        init_guess[4] };
    std::cout << "Tie line" << std::endl;
    bool found_ceq
        = cafe.computeTieLine(temperature, nominalc[0], nominalc[1], &lceq[0]);
    if (found_ceq)
    {
        std::cout << "Found solution: " << lceq[0] << ", " << lceq[1] << ", "
                  << lceq[2] << ", " << lceq[3] << ", " << lceq[4] << "..."
                  << std::endl;
    }

    // use tie line values as initial guesses
    for (int i = 0; i < 4; i++)
        init_guess[i] = lceq[i];

    for (int i = 0; i < nphi + 1; i++)
    {
        double phi = i * dphi;

        double conc[4] = { init_guess[0], init_guess[1], // liquid
            init_guess[2], init_guess[3] }; // solid

        cafe.computePhaseConcentrations(temperature, nominalc, &phi, conc);
        std::cout << "phi=" << phi << ", nominalc=" << nominalc[0] << ", "
                  << nominalc[1] << std::endl;
        std::cout << "KKS solution: cl = (" << conc[0] << "." << conc[1] << ")"
                  << " and cs = (" << conc[2] << "," << conc[3] << ")"
                  << std::endl;
    }
}
