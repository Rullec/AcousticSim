# -------------------enable openmp---------------------
set(USE_OPENMP ON)
if(USE_OPENMP)
    if(WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
        message(STATUS "we use openmp")
    else()
        find_package(OpenMP REQUIRED)
        link_directories(${OpenMP_LIBRARY_DIR})
    endif()
else()
    message(STATUS "we dont use openmp")
endif()

set(CLOTH_SRC
    cloth/BaseCloth.cpp
    cloth/ClothBuilder.cpp
    cloth/BaraffCloth.cpp
    # cloth/BaraffMaterialUnstable.cpp
    cloth/BaraffMaterial.cpp
    cloth/QBendingMaterial.cpp
    cloth/BaseBendingMaterial.cpp
    cloth/DihedralBending.cpp
    cloth/InnerForceCalculator.cpp
    # cloth/ImplicitCloth.cpp   
    # cloth/PBDCloth.cpp
    # cloth/PDCloth.cpp
    # cloth/SemiCloth.cpp
)

set (
    GPU_CLOTH_SRC

    cloth/BaraffClothGpu.cpp
    cloth/BaraffClothGpu.cu
    gpu_utils/CudaMemory.cpp
    gpu_utils/CudaArray.cpp
    gpu_utils/CudaELLMatrix.cpp
    gpu_utils/GPUMatrixOps.cu
    gpu_utils/CudaAsync.cu
    gpu_solvers/PCGSolver.cu
)

if(NOT APPLE)
    list(APPEND CLOTH_SRC ${GPU_CLOTH_SRC})
endif()

# add_executable(cuda_main cloth/BaraffClothGpu.cu)
set(SOFTBODY_SRC
    softbody/SoftBodyImplicit.cpp
    softbody/SoftBody.cpp
    softbody/FourOrderTensor.cpp
    softbody/ThreeOrderTensor.cpp
    softbody/NeoHookeanMaterial.cpp
    softbody/StvkMaterial.cpp
    softbody/MaterialBuilder.cpp
    softbody/BaseMaterial.cpp
    softbody/AcousticLinearElasticity.cpp
    AcousticSoftBody.cpp
    AudioOutput.cpp
    AudioWave.cpp
)
set(
    SIM_SRC
    ${CLOTH_SRC}
    ${SOFTBODY_SRC}
    BaseObject.cpp
    KinematicBody.cpp
    KinematicBodyBuilder.cpp
    Perturb.cpp
    SimObjectBuilder.cpp
    collision/CollisionDetecter.cpp
    collision/BVHCollisionDetecter.cpp
    collision/CollisionInfo.cpp
)

add_library(sim_lib ${SIM_SRC})

# disable warnings of nvcc
# message(STATUS "CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE})
# message(STATUS "CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG})
# message(STATUS "CMAKE_CUDA_FLAGS_DEBUG " ${CMAKE_CUDA_FLAGS_DEBUG})
# message(FATAL_ERROR "CMAKE_CUDA_FLAGS_RELEASE " ${CMAKE_CUDA_FLAGS_RELEASE})

target_compile_options(sim_lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-w>)

pybind11_add_module(softbody_pybind11 ./softbody/softbody_export.cpp)
target_link_libraries(softbody_pybind11 PRIVATE sim_lib utils_lib JSON geo_lib imgui_lib)

add_custom_command(
    TARGET softbody_pybind11 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:softbody_pybind11>
            ${PROJECT_SOURCE_DIR}/dlls)

