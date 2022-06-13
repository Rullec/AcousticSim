cmake_minimum_required(VERSION 3.10)

if(NOT APPLE)
    project(ClothSim VERSION 1.0 LANGUAGES CUDA CXX)
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
    project(ClothSim VERSION 1.0 LANGUAGES CXX)
endif()

add_subdirectory(includes/pybind11)
set(PYBIND11_INCLUDE_PATH "includes/pybind11/include")
include_directories(PYBIND11_INCLUDE_PATH)
include_directories("${PROJECT_SOURCE_DIR}/includes")
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/imgui")


# set version config; the config header will be generated in build dir
configure_file(includes/VersionConfig.h.in VersionConfig.h)
# add build dir into header file

# set language standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_DEBUG_POSTFIX "d")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(configure_movecmd)

if(WIN32)
    link_directories(./lib)
endif()


if(WIN32)
    # add_definitions(-D _USE_MATH_DEFINES)
    # add_definitions(-D _CRT_NO_VA_START_VALIDATION)
    ###### set cpp flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

    ###### set cuda flags
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -G  -Xcompiler=\"-MDd -O0 -Ob0\" ")
    set(CMAKE_CUDA_FLAGS_RELEASE " -Xcompiler=\"-MD -O2 -Ob2\" -DNDEBUG")

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
target_link_libraries(main PUBLIC ${libs} scenes_lib utils_lib geo_lib)
target_include_directories(main PUBLIC ${PROJECT_BINARY_DIR})
if(WIN32)
    move_release_to_src(main)
else()
install(TARGETS main
            RUNTIME
            DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()