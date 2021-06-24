#ifndef included_xlogx

namespace Thermo4PFM
{

double xlogx(const double x);
double xlogx_deriv(const double x);

template <typename DataType>
DataType xlogx_deriv2(const DataType x);
}

#endif
