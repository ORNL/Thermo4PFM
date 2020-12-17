#include "math_utilities.h"
#include <math.h>

namespace Thermo4PFM
{

double Determinant2(double** const m)
{
    return m[0][0] * m[1][1] - m[1][0] * m[0][1];
}

double Determinant3(double** const m)
{
    double d = m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1]
               - m[0][1] * m[1][0] * m[2][2] + m[0][1] * m[1][2] * m[2][0]
               + m[0][2] * m[1][0] * m[2][1] - m[0][2] * m[1][1] * m[2][0];

    return d;
}

//=======================================================================

double Determinant4(double** const m)
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

double Determinant5(double** const m) { return DeterminantN(m, 5); }

//=======================================================================

double DeterminantN(double** mat, const short n)
{
    double* submat[5];
    if (n == 2)
    {
        return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
    }
    else
    {
        double d = 0.;
        for (short c = 0; c < n; c++)
        {
            // loop over rows
            short subi = 0;
            for (short i = 0; i < n; i++)
            {
                if (i == c) continue;

                submat[subi] = &mat[i][1];
                subi++;
            }
            d += (pow(-1, c) * mat[c][0] * DeterminantN(submat, n - 1));
        }
        return d;
    }
}
}
