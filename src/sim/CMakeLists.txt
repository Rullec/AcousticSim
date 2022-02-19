set(
    SIM_SRC
    BaseObject.cpp
    KinematicBody.cpp
    KinematicBodyBuilder.cpp
    Perturb.cpp
    AudioOutput.cpp
    AcousticBody.cpp
    softbody/SoftBody.cpp
    softbody/SoftBodyImplicit.cpp
    softbody/FourOrderTensor.cpp
    softbody/ThreeOrderTensor.cpp
    softbody/NeoHookeanMaterial.cpp
    softbody/BaseMaterial.cpp
    SimObjectBuilder.cpp
)

add_library(sim_lib ${SIM_SRC})