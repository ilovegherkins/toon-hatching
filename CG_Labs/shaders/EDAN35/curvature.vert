#version 410

struct ViewProjTransforms
{
	mat4 view_projection;
	mat4 view_projection_inverse;
};

layout (std140) uniform CameraViewProjTransforms
{
	ViewProjTransforms camera;
};

uniform mat4 normal_model_to_world; 
uniform mat4 vertex_model_to_world;

//uniform sampler2D normals_texture;
//uniform sampler2D tangents_texture;
//uniform sampler2D binormals_texture;


out VS_OUT {
    //vec4 world_pos;
	vec3 normal;
	vec2 texcoord;
	vec3 tangent;
	vec3 binormal;

    //vec3 d_nx;
    //vec3 d_ny;
} vs_out;

//location?
layout (location = 0) in vec3 vertex;
layout (location = 1) in vec3 normal; 
layout (location = 2) in vec3 texcoord;
layout (location = 3) in vec3 tangent;
layout (location = 4) in vec3 binormal;

void main() {
    //vs_out.world_pos = vertex_model_to_world * vec4(vertex, 1.0); // need?
    vs_out.normal = normalize(normal);
    vs_out.binormal = normalize(binormal);
    vs_out.tangent = normalize(tangent);
    vs_out.texcoord = texcoord.xy;

    gl_Position = camera.view_projection * vertex_model_to_world * vec4(vertex, 1.0);
}