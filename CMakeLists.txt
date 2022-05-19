
cmake_minimum_required(VERSION 3.0.0)



if(NOT APPLE)
    project(ClothSim VERSION 0.1.0 LANGUAGES CUDA CXX)
    enable_language(CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    message(STATUS "cuda compiler: " ${CMAKE_CUDA_COMPILER})
    message(STATUS "cuda compiler id: " ${CMAKE_CUDA_COMPILER})
    message(STATUS "cuda compiler version: " ${CMAKE_CUDA_COMPILER_VERSION})

    # cuda CXX version
    if(NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 11)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()
else()
    project(ClothSim VERSION 0.1.0 LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_DEBUG_POSTFIX "d")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(configure_movecmd)
include_directories(./includes)
include_directories(./imgui)
include_directories(./src)

set(DLL_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/dlls)
if(NOT EXISTS ${DLL_OUTPUT_DIR})
    FILE(MAKE_DIRECTORY ${DLL_OUTPUT_DIR})
    message(STATUS "create dll output dir ${DLL_OUTPUT_DIR}")
endif()

if(WIN32)
    link_directories(./lib)
endif()


if(WIN32)
    # add_definitions(-D _USE_MATH_DEFINES)
    # add_definitions(-D _CRT_NO_VA_START_VALIDATION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4061 /wd4244 /wd4246 /wd4305 /wd4267 /wd4711 /wd4710 /wd4514 /wd4477")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /Zi")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
    # use sanitzer
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fsanitize=address")
    # add_definitions("/fsanitize=address")
else()
    # add_definitions(-std=c++17)
    # message(FATAL_ERROR add unix c++17 definitions)
    # ---------------------enable debug symbols---------------------
    add_definitions(-g)
endif()



# ---------------------add address sanitzer---------------------
# if(WIN32)
#     message(STATUS "add address sanitizer")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fsanitize=address")
# else()
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
#     set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
# endif()
find_package(Vulkan REQUIRED)
if(Vulkan_FOUND)
    include_directories(${Vulkan_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "vulkan not found")
endif()

add_subdirectory(./imgui)
add_subdirectory(./src/utils)
add_subdirectory(./src/geometries)
add_subdirectory(./src/sim)
add_subdirectory(./src/cameras)
add_subdirectory(./src/scenes)

set(my_libs ${scenes_lib} ${utils_lib} ${imgui_lib})

add_executable(main ./src/main.cpp)
target_link_libraries(main ${libs} scenes_lib utils_lib geo_lib)

if(WIN32)
    move_release_to_src(main)
else()
install(TARGETS main
            RUNTIME
            DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()