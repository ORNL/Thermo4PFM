#ifndef InterpolationType_H
#define InterpolationType_H

enum class ConcInterpolationType
{
    LINEAR,
    PBG,
    HARMONIC,
    UNDEFINED
};

enum class EnergyInterpolationType
{
    LINEAR,
    PBG,
    HARMONIC,
    UNDEFINED
};

char energyInterpChar(EnergyInterpolationType interp_func_type);
char concInterpChar(ConcInterpolationType interp_func_type);

#endif
