# -------------------enable openmp---------------------
set(USE_OPENMP ON)
if(USE_OPENMP)
    if(WIN32)
        add_definitions(/openmp -DUSE_OPENMP)
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
    cloth/BaraffMaterial.cpp
    cloth/QBendingMaterial.cpp
    # cloth/BaraffMaterialOld.cpp
    # cloth/ImplicitCloth.cpp
    # cloth/PBDCloth.cpp
    # cloth/PDCloth.cpp
    # cloth/SemiCloth.cpp
)

set(SOFTBODY_SRC
    # softbody/SoftBodyImplicit.cpp
    # softbody/SoftBody.cpp
    # softbody/FourOrderTensor.cpp
    # softbody/ThreeOrderTensor.cpp
    # softbody/NeoHookeanMaterial.cpp
    # softbody/StvkMaterial.cpp
    # softbody/MaterialBuilder.cpp
    # softbody/BaseMaterial.cpp
    # AcousticBody.cpp
    # AcousticSoftBody.cpp
)
set(
    SIM_SRC
    ${CLOTH_SRC}
    ${SOFTBODY_SRC}
    BaseObject.cpp
    KinematicBody.cpp
    KinematicBodyBuilder.cpp
    Perturb.cpp
    # AudioOutput.cpp
    SimObjectBuilder.cpp
)

add_library(sim_lib ${SIM_SRC})
# target_link_libraries(sim_lib OpenMP::OpenMP_CXX)