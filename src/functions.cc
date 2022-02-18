#include "functions.h"

#include <math.h>

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
static double (*fun_ptr_arr[3])(
    const double){ linear_interp_func, pbg_interp_func, harmonic_interp_func };

double interp_func(EnergyInterpolationType type, const double phi)
{
    return fun_ptr_arr[static_cast<int>(type)](phi);
}

double interp_func(ConcInterpolationType type, const double phi)
{
    return fun_ptr_arr[static_cast<int>(type)](phi);
}

double pbg_interp_func(const double phi)
{
    double phit = fmax(0., fmin(1., phi));
    return phit * phit * phit * (10. - 15. * phit + 6. * phit * phit);
}

double harmonic_interp_func(const double phi)
{
    double phit = fmax(0., fmin(1., phi));
    return phit * phit * (3. - 2. * phit);
}

double linear_interp_func(const double phi) { return fmax(0., fmin(1., phi)); }

double deriv_pbg_interp_func(const double phi)
{
    double phit = fmax(0., fmin(1., phi));
    return 30. * phit * phit * (1. - phit) * (1. - phit);
}

double deriv_harmonic_interp_func(const double phi)
{
    double phit = fmax(0., fmin(1., phi));
    return 6. * phit * (1. - phit);
}

double deriv_linear_interp_func(const double phi)
{
    if (phi > 0. || phi < 1.)
        return 1.;
    else
        return 0.;
}

double second_deriv_pbg_interp_func(const double phi)
{
    double phit = fmax(0., fmin(1., phi));
    return 60. * phit * (1. - 3. * phit + 2. * phit * phit);
}

double second_deriv_harmonic_interp_func(const double phi)
{
    double phit = fmax(0., fmin(1., phi));
    return 6. * (1. - 2. * phit);
}

double second_deriv_linear_interp_func(const double phi) { return 0.; }
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
