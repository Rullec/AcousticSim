#ifndef PCG_SOLVER_H_
#include "sim/cloth/gpu_utils/Cuda2DArray.h"
#include "sim/cloth/gpu_utils/CudaArray.h"
// class cPCGSolver : public std::enable_shared_from_this<cPCGSolver>
// {
// public:
//     void Solve(const cCudaArray<tCudaVector32i>
//                    &ELL_local_vertex_id_to_global_vertex_id,
//                const cCuda2DArray<tCudaMatrix3f> &A,
//                const cCudaArray<tCudaVector3f> &b,
//                cCudaArray<tCudaVector3f> &x);

// protected:
//     cCudaArray<tCudaVector3f> pcg_rbuf;
//     cCudaArray<tCudaVector3f> pcg_dbuf;
//     cCudaArray<tCudaVector3f> pcg_zbuf;
//     cCudaArray<float> pcg_rMinvr_array;
//     cCudaArray<tCudaMatrix3f> Ainv_precond;
//     cCudaArray<float> pcg_dTAd;
//     std::vector<tCudaVector3f> r_cpu;
// };
#endif