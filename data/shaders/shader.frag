#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(location = 0) in vec4 fragColor;
layout(location = 1) in vec2 fragTexCoord;
layout(location = 2) in vec3 fragCameraPos;
layout(location = 3) in vec3 fragVertexWorldPos;
layout(location = 4) in vec3 fragNormal;

layout(location = 0) out vec4 outColor;

layout(binding = 1) uniform sampler2D texSampler;

bool JudgeIsNan(float val)
{
  return ( val < 0.0 || 0.0 < val || val == 0.0 ) ? false : true;
}

float calc_specular(vec3 light_pos, vec3 pixel_world_pos, vec3 normal, vec3 eye_pos)
{
    // input light L
    vec3 L = normalize(light_pos - pixel_world_pos);
    // output light R
    vec3 R = normalize( 2 * (dot(normal, L)) * normal - L);
    vec3 V = normalize( eye_pos - pixel_world_pos);
    return pow( abs( dot(V, R)), 512);
}
void main() {
    vec3 light_pos = vec3(10.0, 10.0, 10.0);
    vec3 light_color = vec3(1.0, 1.0, 1.0) * 0.5;
    float diffuse_weight = 0.8,
        specular_weight = calc_specular(light_pos, fragVertexWorldPos, fragNormal, fragCameraPos);

    // float sum = diffuse_weight + specular_weight;
    // diffuse_weight /= sum;
    // specular_weight /= sum;
    if(isnan(fragTexCoord.x) == true || isnan( fragTexCoord.y) == true )
    {
        // outColor = vec4(fragColor, 1.0);
        // outColor = vec4(abs(fragNormal), 1.0);
        // outColor = vec4(0.5, 0.5, 1.0, 1.0);
        /*
        phong model:
        I = ka * Ia + kd * Id * (N dot L) + ks * Is * (N dot H)^n
        */
        // normal is nan, only pure color
        if(isnan(fragNormal.x) == true )
        {
            outColor = fragColor;
        }
        else {
            outColor = vec4(specular_weight * light_color + diffuse_weight * vec3(fragColor.x, fragColor.y, fragColor.z), 1.0);
        }

        
    }
    else
    {
        outColor = texture(texSampler, fragTexCoord) + vec4(light_color * specular_weight, 1.0);

        // add blur for remote textures
        float length = length(fragCameraPos - fragVertexWorldPos);
        float begin_blur_dist = 20;
        if(length > begin_blur_dist)
        {   
            outColor += 0.1 * vec4(1, 1, 1 ,1)* log(length - begin_blur_dist + 1);
        }
        
    }
    // outColor[0] = fragColor[3];
    // outColor[1] = fragColor[3];
    // outColor[2] = fragColor[3];
    // outColor[3] = fragColor[3];
    outColor[3] = fragColor[3];
    // outColor[3] = 0.5;
}