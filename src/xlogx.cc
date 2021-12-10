// C2 extension of x(log(x)) function for x<=SMALLX
//
#include <math.h>

#define SMALLX 1.0e-5
#define INVSMALLX 1.e5
#define LOGSMALLX -11.512925464970229
#define SMALLXLOGSMALLX -11.512925464970229e-5

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double xlogx(const double x)
{
    if (x > SMALLX)
    {
        return x * log(x);
    }
    else
    {
        return SMALLXLOGSMALLX + (x - SMALLX) * LOGSMALLX
               + 0.5 * (x * x * INVSMALLX - SMALLX);
    }
}

double xlogx_deriv(const double x)
{
    if (x > SMALLX)
    {
        return log(x) + 1.0;
    }
    else
    {
        return LOGSMALLX + x * INVSMALLX;
    }
}

template <typename DataType>
DataType xlogx_deriv2(const DataType x)
{
    if (x > SMALLX)
    {
        return 1. / x;
    }
    else
    {
        return INVSMALLX;
    }
}

template float xlogx_deriv2<float>(const float);
template double xlogx_deriv2<double>(const double);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
