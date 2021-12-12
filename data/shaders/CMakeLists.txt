
# # add_custom_command(
# #     OUTPUT spv_shaders_linux
# #     COMMAND "sh ./compile.sh"
# #     MAIN_DEPENDENCY shader.vert shader.frag
# #     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
# # )
# # set(
# #     SHADER_SRC
# #     ./shader.frag 
# #     ./shader.vert
# # )
# # add_custom_target(shaders ALL DEPENDS shader.frag shader.vert)

# # add_executable(vulkan_shaders ${SHADER_SRC})
# set(KERNELS
# ${CMAKE_SOURCE_DIR}/src/shaders/shader.vert
# ${CMAKE_SOURCE_DIR}/src/shaders/shader.frag
# )

# set(COMPILED_KERNELS
# ${CMAKE_SOURCE_DIR}/src/shaders/shader.vert.spv
# ${CMAKE_SOURCE_DIR}/src/shaders/shader.frag.spv
# )

# add_executable(vulkan_shaders ${KERNELS} ${COMPILED_KERNELS})

# foreach(KERNEL ${KERNELS})
# add_custom_command(OUTPUT ${KERNEL}.spv
# COMMAND glslc ${KERNEL} -g -o ${KERNEL}.spv
# DEPENDS ${KERNEL}
# COMMENT "Rebuilding ${KERNEL}.spv" )
# message(STATUS "Generating build commands for ${KERNEL}.spv")
# endforeach()