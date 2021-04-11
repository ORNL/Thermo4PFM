
#include "CALPHADFunctions.h"
#include "InterpolationType.h"
#include "functions.h"
#include "xlogx.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <string>

#include <omp.h>

namespace pt = boost::property_tree;

using namespace Thermo4PFM;

int main(int argc, char* argv[])
{
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;

    const double cmin     = 0.1;
    const double cmax     = 0.9;
    const int nCintervals = 10;
    const double deltaC   = (cmax - cmin) / (double)nCintervals;

    std::vector<double> cl(nCintervals + 1);
    std::vector<double> cs(nCintervals + 1);

    {
        // serial loop
        for (int i = 0; i < nCintervals + 1; i++)
        {
            const double c = cmin + i * deltaC;

            // compute concentrations in each phase
            double val = xlogx(c);

            std::cout << "Concentration = " << c << ", xlogx = " << val
                      << std::endl;
        }
    }

// parallel loop
#pragma omp target
#pragma omp parallel for
    for (int i = 0; i < nCintervals + 1; i++)
    {
        const double c = cmin + i * deltaC;

        // compute concentrations in each phase
        double val = xlogx(c);

        double val2 = interp_func(ConcInterpolationType::PBG, 0.5);

        double val3 = CALPHADcomputeFMixBinary(1.2, 1.5, 1.6, 1.7, 0.5);

        //        REQUIRE(conc[0] == Approx(cl[i]).margin(1.e-1));
        //        REQUIRE(conc[1] == Approx(cs[i]).margin(1.e-1));
    }
}
