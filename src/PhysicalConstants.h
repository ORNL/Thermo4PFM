#ifndef included_PhysicalConstants
#define included_PhysicalConstants

namespace Thermo4PFM
{
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
const double gas_constant_R_JpKpmol = 8.314472; // J K-1 mol-1
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}

#endif
