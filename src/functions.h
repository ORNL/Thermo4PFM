#include "InterpolationType.h"

namespace Thermo4PFM
{
double pbg_interp_func(const double phi);
double harmonic_interp_func(const double phi);
double linear_interp_func(const double phi);
double deriv_pbg_interp_func(const double phi);
double deriv_harmonic_interp_func(const double phi);
double deriv_linear_interp_func(const double phi);
double second_deriv_pbg_interp_func(const double phi);
double second_deriv_harmonic_interp_func(const double phi);
double second_deriv_linear_interp_func(const double phi);

double interp_func(EnergyInterpolationType type, const double phi);
double interp_func(ConcInterpolationType type, const double phi);
}
