
#include "Determinant.h"

#include <iostream>
#include <string>

#include <omp.h>

using namespace Thermo4PFM;

#define N 2

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

#ifdef _OPENMP
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;
#endif

    const int nDet = 10;

    double matrix0[2];
    double matrix1[2];

    double* matrix[N] = { &matrix0[0], &matrix1[0] };

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            matrix[i][j] = 2. + rand() / (double)RAND_MAX;

    {
        // serial loop
        for (int i = 0; i < nDet; i++)
        {
            // compute concentrations in each phase
            double val = evalDeterminant<N>(matrix);

            std::cout << "Determinant = " << val << std::endl;
        }
    }

    double* values = new double[nDet];
// parallel loop
// clang-format off
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp target data map(to : matrix0 [0:N]) \
                        map(to : matrix1 [0:N]) \
                        map(from : values [0:nDet])
#endif
// clang-format on
#pragma omp parallel for
    for (int i = 0; i < nDet + 1; i++)
    {
        double* matrix[N] = { &matrix0[0], &matrix1[0] };
        values[i]         = evalDeterminant<N>(matrix);

        //        REQUIRE(conc[0] == Approx(cl[i]).margin(1.e-1));
        //        REQUIRE(conc[1] == Approx(cs[i]).margin(1.e-1));
    }

    // print GPU results
    std::cout << "GPU results:" << std::endl;
    for (int i = 0; i < nDet; i++)
        std::cout << "Determinant = " << values[i] << std::endl;

    delete[] values;
}
