#include <algorithm>
#include <iostream>

double well_func(const double phi)
{
    return 16. * phi * phi * (1. - phi) * (1. - phi);
}

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
