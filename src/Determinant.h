#ifndef Thermo4PFM_included_Determinant
#define Thermo4PFM_included_Determinant

namespace Thermo4PFM
{
template <int N, typename ScalarType>
ScalarType evalDeterminant(ScalarType** const matrix);

template <typename ScalarType>
ScalarType evalDeterminant2(ScalarType** const matrix);

template <typename ScalarType>
ScalarType evalDeterminant3(ScalarType** const matrix);

template <typename ScalarType>
ScalarType evalDeterminant4(ScalarType** const matrix);
}
#endif
