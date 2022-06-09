# add_library(draw_scene_lib

# )

# add_library(sim_scenes_lib 

# )

find_package(GLFW3 REQUIRED)


add_library(
    scenes_lib
    SceneBuilder.cpp
    Scene.cpp
    SimScene.cpp
    SimStateMachine.cpp
    DrawScene.cpp 
    DrawSceneVk.cpp 
    DrawSceneShot.cpp 
    DrawSceneImGUI.cpp
)

target_link_libraries(scenes_lib PRIVATE cam_lib imgui_lib sim_lib ${Vulkan_LIBRARIES} ${GLFW3_LIBRARY})