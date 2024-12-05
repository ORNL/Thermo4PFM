// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "CALPHADFreeEnergyFunctionsBinary.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cerr
            << "ERROR: program needs 2 argument (database + 2 temperatures "
            << std::endl;
        return 1;
    }

    std::string calphad_filename(argv[1]);
    const double temperature_low  = atof(argv[2]);
    const double temperature_high = atof(argv[3]);
    // initial guesses
    double init_guess[2];
    init_guess[0] = atof(argv[4]);
    init_guess[1] = atof(argv[5]);

    {
        Thermo4PFM::EnergyInterpolationType energy_interp_func_type
            = Thermo4PFM::EnergyInterpolationType::PBG;
        Thermo4PFM::ConcInterpolationType conc_interp_func_type
            = Thermo4PFM::ConcInterpolationType::PBG;

        boost::property_tree::ptree calphad_pt;
        assert(calphad_filename.compare(calphad_filename.size() - 4, 4, "json")
               == 0);
        boost::property_tree::read_json(calphad_filename, calphad_pt);

        pt::ptree newton_pt;

        Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(calphad_pt, newton_pt,
            energy_interp_func_type, conc_interp_func_type);

        const int nT = 50;

        std::vector<double> c0s(nT);
        std::vector<double> c0l(nT);

        std::vector<double> fs(nT);
        std::vector<double> fl(nT);

        std::vector<double> d2fs(nT);
        std::vector<double> d2fl(nT);

        double dT = (temperature_high - temperature_low) / 50;

        const double tol = 1.e-8;
        const int maxits = 50;

        // loop over temperature range
        for (int iT = 0; iT < nT; iT++)
        {

            double temperature = temperature_low + iT * dT;
            std::clog << "T = " << temperature << std::endl;
            // compute minimum using Newton iteration to solve for f'(x)=0
            {
                double conc = init_guess[0];
                for (int it = 0; it < maxits; it++)
                {
                    double val = 0.;
                    cafe.computeDerivFreeEnergy(temperature, &conc,
                        Thermo4PFM::PhaseIndex::phaseL, &val);
                    double deriv = 0.;
                    cafe.computeSecondDerivativeFreeEnergy(temperature, &conc,
                        Thermo4PFM::PhaseIndex::phaseL, &deriv);
                    const double delta = val / deriv;
                    conc -= delta;
                    std::clog << conc << std::endl;

                    if (std::abs(delta) < tol)
                    {
                        std::clog << "converged!" << std::endl;
                        d2fl[iT] = deriv;
                        fl[iT]   = cafe.computeFreeEnergy(
                            temperature, &conc, Thermo4PFM::PhaseIndex::phaseL);
                        break;
                    }
                }
                c0l[iT] = conc;
            }

            {
                double conc = init_guess[1];
                for (int it = 0; it < maxits; it++)
                {
                    double val = 0.;
                    cafe.computeDerivFreeEnergy(temperature, &conc,
                        Thermo4PFM::PhaseIndex::phaseA, &val);
                    double deriv = 0.;
                    cafe.computeSecondDerivativeFreeEnergy(temperature, &conc,
                        Thermo4PFM::PhaseIndex::phaseA, &deriv);
                    const double delta = val / deriv;
                    conc -= delta;
                    std::clog << conc << std::endl;

                    if (std::abs(delta) < tol)
                    {
                        std::clog << "converged!" << std::endl;
                        d2fs[iT] = deriv;
                        fs[iT]   = cafe.computeFreeEnergy(
                            temperature, &conc, Thermo4PFM::PhaseIndex::phaseA);
                        break;
                    }
                }
                c0s[iT] = conc;
            }
        }

        {
            std::ofstream os("QuadraticDataLiquid.csv");
            os << "T, c0, f(c0), d2f\n";
            for (int iT = 0; iT < nT; iT++)
            {
                double temperature = temperature_low + iT * dT;
                os << temperature << ", " << c0l[iT] << ", " << fl[iT] << ", "
                   << d2fl[iT] << std::endl;
            }
        }
        {
            std::ofstream os("QuadraticDataSolid.csv");
            os << "T, c0, f(c0), d2f\n";
            for (int iT = 0; iT < nT; iT++)
            {
                double temperature = temperature_low + iT * dT;
                os << temperature << ", " << c0s[iT] << ", " << fs[iT] << ", "
                   << d2fs[iT] << std::endl;
            }
        }
    }

    return 0;
}
