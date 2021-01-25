#include "Determinant.h"

#include <math.h>

namespace Thermo4PFM
{

template <>
double evalDeterminant<2>(double** const m)
{
    return m[0][0] * m[1][1] - m[1][0] * m[0][1];
}

template <>
double evalDeterminant<3>(double** const m)
{
    double d = m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1]
               - m[0][1] * m[1][0] * m[2][2] + m[0][1] * m[1][2] * m[2][0]
               + m[0][2] * m[1][0] * m[2][1] - m[0][2] * m[1][1] * m[2][0];

    return d;
}

//=======================================================================

template <>
double evalDeterminant<4>(double** const m)
{
    double d
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

template <int N>
double evalDeterminant(double** mat)
{
    double* submat[N - 1];

    double d = 0.;
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
        d += (minus1pow(c) * mat[c][0] * evalDeterminant<N - 1>(submat));
    }
    return d;
}

template double evalDeterminant<5>(double**);
}
