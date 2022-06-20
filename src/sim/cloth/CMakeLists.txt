set(CLOTH_SRC
BaseCloth.cpp
ClothBuilder.cpp
BaraffCloth.cpp
# BaraffMaterialUnstable.cpp
BaraffMaterial.cpp
QBendingMaterial.cpp
BaseBendingMaterial.cpp
DihedralBending.cpp
InnerForceCalculator.cpp
# ImplicitCloth.cpp   
# PBDCloth.cpp
# PDCloth.cpp
# SemiCloth.cpp
)

set (
GPU_CLOTH_SRC

BaraffClothGpu.cpp
BaraffClothGpu.cu
gpu_utils/CudaMemory.cpp
gpu_utils/CudaArray.cpp
gpu_utils/CudaELLMatrix.cpp
gpu_utils/GPUMatrixOps.cu
gpu_utils/CudaAsync.cu
gpu_solvers/PCGSolver.cu
)

# do not compile cuda code on APPLE
if(NOT APPLE)
    list(APPEND CLOTH_SRC ${GPU_CLOTH_SRC})
endif()

add_library(cloth_lib ${CLOTH_SRC})

target_compile_options(cloth_lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-w>)
