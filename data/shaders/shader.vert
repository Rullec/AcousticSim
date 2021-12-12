#version 450
#extension GL_ARB_separate_shader_objects : enable


layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec4 inColor;
layout(location = 2) in vec3 inNormal;
layout(location = 3) in vec2 inTexCoord;

layout(location = 0) out vec4 fragColor;
layout(location = 1) out vec2 fragTexCoord;
layout(location = 2) out vec3 fragCameraPos;
layout(location = 3) out vec3 fragVertexWorldPos;
layout(location = 4) out vec3 fragNormal;

layout(binding = 0) uniform UniformBufferObject {
    vec4 camera_pos;
    mat4 model;
    mat4 view;
    mat4 proj;
} ubo;

void main() {
    gl_Position = ubo.proj * ubo.view * ubo.model * vec4(inPosition, 1.0);
    fragColor = inColor;
    fragTexCoord = inTexCoord;

    // output camera pos and vertex world pos
    fragCameraPos = vec3(ubo.camera_pos.xyz);
    fragVertexWorldPos = vec3((ubo.model * vec4(inPosition, 1.0)).xyz);
    fragNormal = inNormal;
}