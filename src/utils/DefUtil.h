#include <memory>

// Predefined
#define SIM_DECLARE_PTR(ClassObject)                                           \
    using ClassObject##Ptr = std::shared_ptr<ClassObject>;

#define SIM_DECLARE_CLASS_AND_PTR(ClassObject)                                 \
    class ClassObject;                                                         \
    SIM_DECLARE_PTR(ClassObject);

#define SIM_DECLARE_STRUCT_AND_PTR(ClassObject)                                \
    struct ClassObject;                                                        \
    SIM_DECLARE_PTR(ClassObject);

// posx3, colorx4, normalx3, uvx2
#define RENDERING_SIZE_PER_VERTICE (12)


#define SIM_MAX(a, b) ( (a > b) ? a : b)
#define SIM_MIN(a, b) ( (a > b) ? b : a)
#define SIM_SWAP(a, b) {auto c = a; a = b; b = c;}
/*************************************************************************
***************************       OpenMP      ****************************
*************************************************************************/

#ifdef EnableOMP
#include <omp.h>
#define OMP_NUM_THREADS ((omp_get_num_procs()-1)>1?(omp_get_num_procs()-1):1)
#ifdef _WIN32
#define OMP_BARRIER __pragma(omp barrier)
#define OMP_PARALLEL __pragma(omp parallel num_threads(OMP_NUM_THREADS))
#define OMP_PARALLEL_FOR __pragma(omp parallel for num_threads(OMP_NUM_THREADS))
#define OMP_PARALLEL_FOR_SUM_REDUCTION(sum) __pragma(omp parallel for num_threads(OMP_NUM_THREADS) reduction(+: sum))
#else
#define OMP_BARRIER _Pragma("omp barrier")
#define OMP_PARALLEL _Pragma("omp parallel num_threads(OMP_NUM_THREADS)")
#define OMP_PARALLEL_FOR _Pragma("omp parallel for num_threads(OMP_NUM_THREADS)")
#define OMP_PARALLEL_FOR_SUM_REDUCTION(sum) _Pragma("omp parallel for num_threads(OMP_NUM_THREADS) reduction(+: sum)")
#endif
#else
#define OMP_NUM_THREADS 1
#define OMP_BARRIER 
#define OMP_PARALLEL 
#define OMP_PARALLEL_FOR 
#define OMP_PARALLEL_FOR_SUM_REDUCTION(sum)
#endif
