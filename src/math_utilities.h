#ifndef included_math_utilities
#define included_math_utilities

namespace Thermo4PFM
{

double Determinant2(double** const m);
double Determinant3(double** const m);
double Determinant4(double** const m);
double Determinant5(double** const m);
double DeterminantN(double** const m, const short n);
}

#endif
