#include "InterpolationType.h"

namespace Thermo4PFM
{

char concInterpChar(ConcInterpolationType interp_func_type)
{
    switch (interp_func_type)
    {
        case ConcInterpolationType::LINEAR:
            return 'l';
        case ConcInterpolationType::PBG:
            return 'p';
        case ConcInterpolationType::HARMONIC:
            return 'h';
        default:
            return '0';
    }
}

char energyInterpChar(EnergyInterpolationType interp_func_type)
{
    switch (interp_func_type)
    {
        case EnergyInterpolationType::LINEAR:
            return 'l';
        case EnergyInterpolationType::PBG:
            return 'p';
        case EnergyInterpolationType::HARMONIC:
            return 'h';
        default:
            return '0';
    }
}
}
