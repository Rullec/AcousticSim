set(
    SIM_SRC
    BaseObject.cpp
    KinematicBody.cpp
    KinematicBodyBuilder.cpp
    Perturb.cpp
)

add_library(sim_lib ${SIM_SRC})