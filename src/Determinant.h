#ifndef included_Determinant
#define included_Determinant

namespace Thermo4PFM
{

class Determinant
{
public:
    Determinant(){};

    template <int N>
    static double evaluate(double** const matrix);
};
}
#endif
