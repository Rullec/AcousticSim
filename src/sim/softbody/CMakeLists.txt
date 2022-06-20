set(SOFTBODY_SRC
    SoftBodyImplicit.cpp
    SoftBody.cpp
    FourOrderTensor.cpp
    ThreeOrderTensor.cpp
    NeoHookeanMaterial.cpp
    StvkMaterial.cpp
    MaterialBuilder.cpp
    BaseMaterial.cpp
    AcousticLinearElasticity.cpp
)

add_library(softbody_lib ${SOFTBODY_SRC})