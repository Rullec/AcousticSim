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

add_subdirectory(acoustic)
add_subdirectory(cloth)
add_subdirectory(collision)
add_subdirectory(kinematic)
add_subdirectory(softbody)

set(
    SIM_BASE_SRC    
    BaseObject.cpp
    Perturb.cpp
    SimObjectBuilder.cpp
    ViscoMassSpring.cpp
)

add_library(sim_lib ${SIM_BASE_SRC})

target_link_libraries(sim_lib PRIVATE 
    acoustic_lib
    cloth_lib
    collision_lib
    kinematic_lib
    softbody_lib
)