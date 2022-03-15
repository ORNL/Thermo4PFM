#define CATCH_CONFIG_MAIN

#include "CALPHADConcSolverTernary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "CALPHADFunctions.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "PhysicalConstants.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD ternary solver", "[ternary solver]")
{
    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphadMoNbTa.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    double temperature = 2923.;

    CalphadDataType LmixABPhaseL[4][2];
    CalphadDataType LmixABPhaseA[4][2];

    CalphadDataType LmixACPhaseL[4][2];
    CalphadDataType LmixACPhaseA[4][2];

    CalphadDataType LmixBCPhaseL[4][2];
    CalphadDataType LmixBCPhaseA[4][2];

    CalphadDataType LmixABCPhaseL[3][2];
    CalphadDataType LmixABCPhaseA[3][2];

    {
        std::string dbnamemixL("LmixABCPhaseL");
        if (calphad_db.get_child_optional(dbnamemixL))
        {
            pt::ptree& Lmix0_db = calphad_db.get_child(dbnamemixL);
            Thermo4PFM::readLmixTernaryParameters(Lmix0_db, LmixABCPhaseL);
        }
        else
        {
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 2; i++)
                {
                    LmixABCPhaseL[j][i] = 0.;
                }
        }
    }
    {
        std::string dbnamemixL("LmixABCPhaseA");
        if (calphad_db.get_child_optional(dbnamemixL))
        {
            pt::ptree& Lmix0_db = calphad_db.get_child(dbnamemixL);
            Thermo4PFM::readLmixTernaryParameters(Lmix0_db, LmixABCPhaseA);
        }
        else
        {
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 2; i++)
                {
                    LmixABCPhaseA[j][i] = 0.;
                }
        }
    }

    {
        std::string dbnamemixL("LmixABPhaseL");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixL);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixABPhaseL);
    }
    {
        std::string dbnamemixA("LmixABPhaseA");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixA);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixABPhaseA);
    }
    {
        std::string dbnamemixL("LmixACPhaseL");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixL);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixACPhaseL);
    }
    {
        std::string dbnamemixA("LmixACPhaseA");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixA);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixACPhaseA);
    }
    {
        std::string dbnamemixL("LmixBCPhaseL");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixL);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixBCPhaseL);
    }
    {
        std::string dbnamemixA("LmixBCPhaseA");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixA);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixBCPhaseA);
    }

    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy g_species_phaseL[3];
    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy g_species_phaseA[3];
    {
        std::string dbnameL("PhaseL");
        std::string dbnameA("PhaseA");

        pt::ptree& speciesA_db = calphad_db.get_child("SpeciesA");
        g_species_phaseL[0].initialize("L0", speciesA_db.get_child(dbnameL));
        g_species_phaseA[0].initialize("A0", speciesA_db.get_child(dbnameA));

        pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
        g_species_phaseL[1].initialize("L1", speciesB_db.get_child(dbnameL));
        g_species_phaseA[1].initialize("A1", speciesB_db.get_child(dbnameA));

        pt::ptree& speciesC_db = calphad_db.get_child("SpeciesC");
        g_species_phaseL[2].initialize("L2", speciesC_db.get_child(dbnameL));
        g_species_phaseA[2].initialize("A2", speciesC_db.get_child(dbnameA));
    }

    CalphadDataType fA[2];
    fA[0] = g_species_phaseL[0].fenergy(temperature);
    fA[1] = g_species_phaseA[0].fenergy(temperature);
    // std::cout<<"fA[0]="<<fA[0]<<", fA[1]="<<fA[1]<<std::endl;

    CalphadDataType fB[2];
    fB[0] = g_species_phaseL[1].fenergy(temperature);
    fB[1] = g_species_phaseA[1].fenergy(temperature);
    // std::cout<<"fB[0]="<<fB[0]<<", fB[1]="<<fB[1]<<std::endl;

    CalphadDataType fC[2];
    fC[0] = g_species_phaseL[2].fenergy(temperature);
    fC[1] = g_species_phaseA[2].fenergy(temperature);

    CalphadDataType L_AB_L[4];
    for (int i = 0; i < 4; i++)
        L_AB_L[i] = LmixABPhaseL[i][0] + temperature * LmixABPhaseL[i][1];
    CalphadDataType L_AB_S[4];
    for (int i = 0; i < 4; i++)
        L_AB_S[i] = LmixABPhaseA[i][0] + temperature * LmixABPhaseA[i][1];

    CalphadDataType L_AC_L[4];
    for (int i = 0; i < 4; i++)
        L_AC_L[i] = LmixACPhaseL[i][0] + temperature * LmixACPhaseL[i][1];

    CalphadDataType L_AC_S[4];
    for (int i = 0; i < 4; i++)
        L_AC_S[i] = LmixACPhaseA[i][0] + temperature * LmixACPhaseA[i][1];

    CalphadDataType L_BC_L[4];
    for (int i = 0; i < 4; i++)
        L_BC_L[i] = LmixBCPhaseL[i][0] + temperature * LmixBCPhaseL[i][1];

    CalphadDataType L_BC_S[4];
    for (int i = 0; i < 4; i++)
        L_BC_S[i] = LmixBCPhaseA[i][0] + temperature * LmixBCPhaseA[i][1];

    CalphadDataType L_ABC_L[3];
    for (int i = 0; i < 3; i++)
        L_ABC_L[i] = LmixABCPhaseL[i][0] + temperature * LmixABCPhaseL[i][1];

    CalphadDataType L_ABC_S[3];
    for (int i = 0; i < 3; i++)
        L_ABC_S[i] = LmixABCPhaseA[i][0] + temperature * LmixABCPhaseA[i][1];

    const double RTinv = 1.0 / (GASCONSTANT_R_JPKPMOL * temperature);

    double sol[4] = { 0.33, 0.38, 0.32, 0.33 };
    double hphi   = 0.5;
    double c0     = 0.33;
    double c1     = 0.33;

    Thermo4PFM::CALPHADConcSolverTernary solver;
    solver.setup(c0, c1, hphi, RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S, L_AC_S,
        L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);

    int nits = solver.ComputeConcentration(sol, 1.e-8, 50);

    std::cout << "Solution = " << sol[0] << "," << sol[1] << "," << sol[2]
              << "," << sol[3] << std::endl;
    double ref_sol[4] = { 0.339215, 0.356739, 0.320785, 0.303261 };

    CHECK(sol[0] == Approx(ref_sol[0]).margin(1.e-6));
    CHECK(sol[1] == Approx(ref_sol[1]).margin(1.e-6));
    CHECK(sol[2] == Approx(ref_sol[2]).margin(1.e-6));
    CHECK(sol[3] == Approx(ref_sol[3]).margin(1.e-6));

    std::cout << "nits=" << nits << std::endl;
}
