
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "InterpolationType.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

#include <omp.h>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

#ifdef _OPENMP
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;
#endif

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    const int nintervals = 10;
    double deviation     = 0.1 / (double)nintervals;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphadAuNi.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    double temperature = 1450.;
    std::cout << "Temperature = " << temperature << std::endl;

    // initial guesses
    const double init_guess[2] = { 0.5, 0.5 };

    std::vector<double> cl(nintervals);
    std::vector<double> cs(nintervals);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary* cafe
        = new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(calphad_db,
            newton_db, energy_interp_func_type, conc_interp_func_type);

    short* nits = new short[nintervals];

    {
        // serial loop
        for (int i = 0; i < nintervals; i++)
        {
            double phi = 0.5 + i * deviation;
            double c0  = 0.3;

            // solution
            double conc[2] = { init_guess[0], init_guess[1] };

            // compute concentrations in each phase
            nits[i] = cafe->computePhaseConcentrations(
                temperature, &c0, &phi, conc);

            std::cout << "Number of Newton iterations: " << nits[i]
                      << std::endl;
            std::cout << "Concentrations: cl = " << conc[0]
                      << " and cs = " << conc[1] << "..." << std::endl;

            cl[i] = conc[0];
            cs[i] = conc[1];
        }
    }

#ifdef HAVE_OPENMP_OFFLOAD
    double* xdev = new double[2 * nintervals];
    for (int i = 0; i < 2 * nintervals; i++)
    {
        xdev[i] = -1;
    }

// parallel loop
// clang-format off
#pragma omp target map(to : cafe [0:1]) \
                   map(to : init_guess[:2]) \
                   map(from : xdev) \
                   map(from : nits[:nintervals])
// clang-format on
#pragma omp parallel for
    for (int i = 0; i < nintervals + 1; i++)
    {
        double phi = 0.5 + i * deviation;
        double c0  = 0.3;

        xdev[2 * i]     = init_guess[0];
        xdev[2 * i + 1] = init_guess[1];

        // compute concentrations in each phase
        nits[i] = cafe->computePhaseConcentrations(
            temperature, &c0, &phi, &xdev[2 * i]);
    }

    for (int i = 0; i < nintervals; i++)
    {
        std::cout << "Number of Newton iterations: " << nits[i] << std::endl;
        std::cout << "Device: x=" << xdev[2 * i] << "," << xdev[2 * i + 1]
                  << std::endl;
    }
    delete[] xdev;
#endif
    delete[] nits;
}
