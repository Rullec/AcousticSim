
cmake_minimum_required(VERSION 3.0.0)
project(ClothSim VERSION 0.1.0)

set(CMAKE_DEBUG_POSTFIX "d")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(configure_movecmd)
include_directories(./includes)
include_directories(./imgui)

set(DLL_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/dlls)
if(NOT EXISTS ${DLL_OUTPUT_DIR})
    FILE(MAKE_DIRECTORY ${DLL_OUTPUT_DIR})
    message(STATUS "create dll output dir ${DLL_OUTPUT_DIR}")
endif()

if(WIN32)
    link_directories(./lib)
endif()

find_package(Vulkan REQUIRED)
if(Vulkan_FOUND)
    include_directories(${Vulkan_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "vulkan not found")
endif()


find_package(GLFW3 REQUIRED)
include_directories(./src)
find_package(OpenMP REQUIRED)
link_directories(${OpenMP_LIBRARY_DIR})

if(WIN32)
    add_definitions(-D _USE_MATH_DEFINES)
    add_definitions(-D _CRT_NO_VA_START_VALIDATION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /bigobj")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4061 /wd4244 /wd4246 /wd4305 /wd4267 /wd4711 /wd4710 /wd4514 /wd4477")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /Zi")
    add_definitions(/arch:AVX2)
    # use sanitzer
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fsanitize=address")
    # add_definitions("/fsanitize=address")
else()
    add_definitions(-std=c++17)
    # message(FATAL_ERROR add unix c++17 definitions)
    # ---------------------enable debug symbols---------------------
    add_definitions(-g)
endif()


# -------------------enable openmp---------------------
set(USE_OPENMP ON)
if(USE_OPENMP)
    add_definitions(/openmp -DUSE_OPENMP)
    message(STATUS "we use openmp")
else()
    message(STATUS "we dont use openmp")
endif()

# -------------------enable optix-----------------
if(USE_OPTIX)
    include(./cmake/configure_optix.cmake)
endif()

# ---------------------add address sanitzer---------------------
# if(WIN32)
#     message(STATUS "add address sanitizer")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fsanitize=address")
# else()
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
#     set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
# endif()


add_subdirectory(./imgui)
add_subdirectory(./src/utils)
add_subdirectory(./src/geometries)
add_subdirectory(./src/sim)
add_subdirectory(./src/cameras)
add_subdirectory(./src/scenes)

set(my_libs ${scenes_lib} ${utils_lib} ${imgui_lib})
set(libs ${Vulkan_LIBRARIES} ${GLFW3_LIBRARY} ${my_libs})

add_executable(main ./src/main.cpp ./src/CGSolver.cpp)
target_link_libraries(main ${libs} scenes_lib utils_lib geo_lib OpenMP::OpenMP_CXX)

if(WIN32)
    install(TARGETS main
            RUNTIME
            DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    )
else()
    move_release_to_src(main)
endif()