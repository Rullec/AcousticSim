
add_library(
    geo_lib
    Primitives.cpp
    Triangulator.cpp
    Raycaster.cpp
    BezierCurve.cpp
    Tetrahedron.cpp
    TriangulatorDelaunay.cpp
    ObjectBVH.cpp
    AABB.cpp
    Arrow.cpp
)

target_link_libraries(geo_lib PRIVATE utils_lib)