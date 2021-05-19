// C2 extension of x(log(x)) function for x<=smallx
//
#include <math.h>

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
static const double smallx            = 1.0e-5;
static const double inv_smallx        = 1. / smallx;
static const double log_smallx        = log(smallx);
static const double smallx_log_smallx = smallx * log_smallx;
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double xlogx(const double x)
{
    if (x > smallx)
    {
        return x * log(x);
    }
    else
    {
        return smallx_log_smallx + (x - smallx) * log_smallx
               + 0.5 * (x * x * inv_smallx - smallx);
    }
}

double xlogx_deriv(const double x)
{
    if (x > smallx)
    {
        return log(x) + 1.0;
    }
    else
    {
        return log_smallx + x * inv_smallx;
    }
}

double xlogx_deriv2(const double x)
{
    if (x > smallx)
    {
        return 1. / x;
    }
    else
    {
        return inv_smallx;
    }
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
