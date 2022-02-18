add_library(draw_scene_lib
DrawScene.cpp 
DrawSceneVk.cpp 
DrawSceneShot.cpp 
DrawSceneImGUI.cpp
)

add_library(sim_scenes_lib 
Scene.cpp
SimScene.cpp
SimStateMachine.cpp
)


add_library(
    scenes_lib
    SceneBuilder.cpp
)

target_link_libraries(sim_scenes_lib sim_lib geo_lib utils_lib cam_lib)
target_link_libraries(draw_scene_lib ${Vulkan_LIBRARIES} imgui_lib cam_lib)
target_link_libraries(scenes_lib ${GLFW3_LIBRARY} draw_scene_lib sim_scenes_lib cam_lib)