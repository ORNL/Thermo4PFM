#include <algorithm>
#include <iostream>

namespace Thermo4PFM
{

double well_func(const double phi)
{
    return 16. * phi * phi * (1. - phi) * (1. - phi);
}

double deriv_well_func(const double phi)
{
    return 32. * phi * (1. - phi) * (1. - 2. * phi);
}

double second_deriv_well_func(const double phi)
{
    return 32. * (1. + 6. * phi * (phi - 1.));
}
}
