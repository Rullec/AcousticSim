#include "CudaMatrix.h"
// #include "CudaVector.h"
#include <math.h>
class cCudaMath
{
public:
    // for unrecognized value, it cannot match the std::isnan, and we return
    // false
    template <typename Type> static bool IsNan(const Type &val)
    {
        return false;
    };

    template <> static bool IsNan(const float &val) { return isnan(val); };
    template <> static bool IsNan(const double &val)
    {
        return isnan(val);
    };

    template <int N, int M>
    SIM_CUDA_CALLABLE_INLINE static bool
    IsNan(const tCudaMatrix<float, N, M> &a)
    {
        for (int i = 0; i < a.mElements; i++)
        {
            if (isnan(a.mData[i]))
                return true;
        }
        return false;
    }
};