#include "CALPHADSpeciesPhaseGibbsEnergy.h"

#include <cassert>
#include <iostream>
#include <math.h>
#include <vector>

#ifndef CALPHAD4PFM_ERROR

#define CALPHAD4PFM_ERROR(X)                                                   \
    do                                                                         \
    {                                                                          \
        std::cerr << "ERROR in file " << __FILE__ << " at line " << __LINE__   \
                  << std::endl;                                                \
        std::cerr << "Error Message: " << X << std::endl;                      \
        abort();                                                               \
    } while (0)

#endif

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

void read_optional(pt::ptree& db, const std::string key,
    std::vector<double>& coeffs, const unsigned int nintervals)
{
    if (db.get_child_optional(key))
    {
        for (pt::ptree::value_type& v : db.get_child(key))
        {
            coeffs.push_back(v.second.get_value<double>());
        }
        assert(nintervals == coeffs.size());
    }
    else
    {
        for (unsigned i = 0; i < nintervals; i++)
            coeffs.push_back(0.);
    }
}

void CALPHADSpeciesPhaseGibbsEnergy::initialize(
    const std::string& name, pt::ptree& db)
{
    name_ = new char[name.length() + 1];
    strcpy(name_, name.c_str());

    std::vector<double> tmp;
    for (pt::ptree::value_type& tc : db.get_child("Tc"))
    {
        tmp.push_back(tc.second.get_value<double>());
    }

    size_t ntc = tmp.size();
    assert(ntc > 1);
    assert(ntc < MAXNINTERVALS + 2);

    for (size_t i = 0; i < ntc; i++)
        tc_[i] = tmp[i];
    nintervals_ = ntc - 1;

    std::vector<double> a;
    for (pt::ptree::value_type& v : db.get_child("a"))
    {
        a.push_back(v.second.get_value<double>());
    }
    const size_t nintervals = a.size();
    assert(nintervals == ntc - 1);

    std::vector<double> b;
    for (pt::ptree::value_type& v : db.get_child("b"))
    {
        b.push_back(v.second.get_value<double>());
    }
    assert(nintervals == b.size());

    std::vector<double> c;
    for (pt::ptree::value_type& v : db.get_child("c"))
    {
        c.push_back(v.second.get_value<double>());
    }
    assert(nintervals == c.size());

    std::vector<double> d2;
    for (pt::ptree::value_type& v : db.get_child("d2"))
    {
        d2.push_back(v.second.get_value<double>());
    }
    assert(nintervals == d2.size());

    std::vector<double> d3;
    read_optional(db, "d3", d3, nintervals);

    std::vector<double> d4;
    read_optional(db, "d4", d4, nintervals);

    std::vector<double> d7;
    read_optional(db, "d7", d7, nintervals);

    std::vector<double> dm1;
    read_optional(db, "dm1", dm1, nintervals);

    std::vector<double> dm9;
    read_optional(db, "dm9", dm9, nintervals);

    if (db.get_child_optional("d5"))
    {
        CALPHAD4PFM_ERROR(
            "CALPHADSpeciesPhaseGibbsEnergy: T**5 not implemented!!!");
    }
    if (db.get_child_optional("d6"))
    {
        CALPHAD4PFM_ERROR(
            "CALPHADSpeciesPhaseGibbsEnergy: T**6 not implemented!!!");
    }
    if (db.get_child_optional("dm2"))
    {
        CALPHAD4PFM_ERROR(
            "CALPHADSpeciesPhaseGibbsEnergy: T**-2 not implemented!!!");
    }

    double* pa   = a.data();
    double* pb   = b.data();
    double* pc   = c.data();
    double* pd2  = d2.data();
    double* pd3  = d3.data();
    double* pd4  = d4.data();
    double* pd7  = d7.data();
    double* pdm1 = dm1.data();
    double* pdm9 = dm9.data();

    for (unsigned i = 0; i < nintervals; i++)
    {
        expansion_[i].init(pa[i], pb[i], pc[i], pd2[i], pd3[i], pd4[i], pd7[i],
            pdm1[i], pdm9[i]);
    }
}

/////////////////////////////////////////////////////////////////////
// Free energy function
// parameters are in J/mol
// returned values are in J/mol
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double CALPHADSpeciesPhaseGibbsEnergy::fenergy(
    const double T) // expect T in Kelvin
{
    // assert(nintervals_ > 0);

    for (int i = 0; i < nintervals_; i++)
        if (T >= tc_[i] && T < tc_[i + 1])
        {
            return expansion_[i].value(T);
        }

    return NAN;
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

void CALPHADSpeciesPhaseGibbsEnergy::plotFofT(
    std::ostream& os, const double T0, const double T1)
{
    const double dT = 10.;
    const int npts  = (int)std::trunc((T1 - T0) / dT);
    std::string name(name_);
    os << "# fenergy(J/mol) vs. T(K) for species " << name << std::endl;
    for (int i = 0; i < npts; i++)
    {
        double testT = T0 + dT * i;
        os << testT << '\t';
        os << fenergy(testT) << std::endl;
    }
    os << std::endl;
}
}
