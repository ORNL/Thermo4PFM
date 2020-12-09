#include <algorithm>
#include <iostream>

namespace Thermo4PFM
{

double interp_func(const double phi, const char type)
{
    switch (type)
    {
        case 'q':
        {
            double phit = std::max(0., phi);
            return phit * phit;
        }
        case 'p':
        {
            double phit = std::max(0., std::min(1., phi));
            return phit * phit * phit * (10. - 15. * phit + 6. * phit * phit);
        }
        case 'h':
        {
            double phit = std::max(0., std::min(1., phi));
            return phit * phit * (3. - 2. * phit);
        }
        case 'l':
        {
            double phit = std::max(0., std::min(1., phi));
            return phit;
        }
        default:
        {
            std::cerr << "Unknown interpolation type" << std::endl;
            abort();
            return 0.;
        }
    }
}

double deriv_interp_func(const double phi, const char type)
{
    switch (type)
    {
        case 'q':
        {
            double phit = std::max(0., phi);
            return 2. * phit;
        }
        case 'p':
        {
            double phit = std::max(0., std::min(1., phi));
            return 30. * phit * phit * (1. - phit) * (1. - phit);
        }
        case 'h':
        {
            double phit = std::max(0., std::min(1., phi));
            return 6. * phit * (1. - phit);
        }
        case 'l':
        {
            if (phi > 0. || phi < 1.)
                return 1.;
            else
                return 0.;
        }
        default:
        {
            std::cerr << "Unknown interpolation type" << std::endl;
            abort();
            return 0.;
        }
    }
}
}
