# add_library(draw_scene_lib

# )

# add_library(sim_scenes_lib 

# )


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

target_link_libraries(scenes_lib cam_lib imgui_lib sim_lib ${Vulkan_LIBRARIES} ${GLFW3_LIBRARY})