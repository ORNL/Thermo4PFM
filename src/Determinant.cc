#include "Determinant.h"

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

template <typename ScalarType>
ScalarType evalDeterminant2(ScalarType** const m)
{
    return m[0][0] * m[1][1] - m[1][0] * m[0][1];
}

template <typename ScalarType>
ScalarType evalDeterminant3(ScalarType** const m)
{
    ScalarType d = m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1]
                   - m[0][1] * m[1][0] * m[2][2] + m[0][1] * m[1][2] * m[2][0]
                   + m[0][2] * m[1][0] * m[2][1] - m[0][2] * m[1][1] * m[2][0];

    return d;
}

//=======================================================================

template <typename ScalarType>
ScalarType evalDeterminant4(ScalarType** const m)
{
    ScalarType d
        = m[0][0]
              * (m[1][1] * m[2][2] * m[3][3] - m[1][1] * m[2][3] * m[3][2]
                    - m[1][2] * m[2][1] * m[3][3] + m[1][2] * m[2][3] * m[3][1]
                    + m[1][3] * m[2][1] * m[3][2] - m[1][3] * m[2][2] * m[3][1])

          - m[0][1]
                * (m[1][0] * m[2][2] * m[3][3] - m[1][0] * m[2][3] * m[3][2]
                      - m[1][2] * m[2][0] * m[3][3]
                      + m[1][2] * m[2][3] * m[3][0]
                      + m[1][3] * m[2][0] * m[3][2]
                      - m[1][3] * m[2][2] * m[3][0])

          + m[0][2]
                * (m[1][0] * m[2][1] * m[3][3] - m[1][0] * m[2][3] * m[3][1]
                      - m[1][1] * m[2][0] * m[3][3]
                      + m[1][1] * m[2][3] * m[3][0]
                      + m[1][3] * m[2][0] * m[3][1]
                      - m[1][3] * m[2][1] * m[3][0])

          - m[0][3]
                * (m[1][0] * m[2][1] * m[3][2] - m[1][0] * m[2][2] * m[3][1]
                      - m[1][1] * m[2][0] * m[3][2]
                      + m[1][1] * m[2][2] * m[3][0]
                      + m[1][2] * m[2][0] * m[3][1]
                      - m[1][2] * m[2][1] * m[3][0]);

    return d;
}

//=======================================================================

// return pow(-1,i)
int minus1pow(const short i) { return i & 1 ? -1 : 1; }

template <>
double evalDeterminant<2, double>(double** m)
{
    return evalDeterminant2(m);
}

template <>
float evalDeterminant<2, float>(float** m)
{
    return evalDeterminant2(m);
}

template <>
double evalDeterminant<3, double>(double** m)
{
    return evalDeterminant3(m);
}

template <>
float evalDeterminant<3, float>(float** m)
{
    return evalDeterminant3(m);
}

template <>
double evalDeterminant<4, double>(double** m)
{
    return evalDeterminant4(m);
}

template <>
float evalDeterminant<4, float>(float** m)
{
    return evalDeterminant4(m);
}

template <int N, typename ScalarType>
ScalarType evalDeterminant(ScalarType** mat)
{
    ScalarType* submat[N - 1];

    ScalarType d = 0.;
    for (short c = 0; c < N; c++)
    {
        // loop over rows
        short subi = 0;
        for (short i = 0; i < N; i++)
        {
            if (i == c) continue;

            submat[subi] = &mat[i][1];
            subi++;
        }
        d += (minus1pow(c) * mat[c][0]
              * evalDeterminant<N - 1, ScalarType>(submat));
    }
    return d;
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

template double evalDeterminant<5, double>(double**);
template float evalDeterminant<5, float>(float**);
}
