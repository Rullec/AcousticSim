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