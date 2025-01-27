#version 410

uniform bool has_diffuse_texture;
uniform bool has_specular_texture;
uniform bool has_normals_texture;
uniform bool has_opacity_texture;
uniform sampler2D diffuse_texture;
uniform sampler2D specular_texture;
uniform sampler2D normals_texture;
uniform sampler2D opacity_texture;
uniform mat4 normal_model_to_world; 
uniform mat4 vertex_model_to_world;
uniform float hatch_spacing;

uniform vec3 camera_position;
in VS_OUT {
	vec4 world_pos;
	vec3 normal;
	vec2 texcoord;
	vec3 tangent;
	vec3 binormal;
} fs_in;

layout (location = 0) out vec4 geometry_diffuse;
layout (location = 1) out vec4 geometry_specular;
layout (location = 2) out vec4 geometry_normal;
layout (location = 3) out vec4 geometry_tangent;
layout (location = 4) out vec4 geometry_binormal; 

void main()
{	
	if (has_opacity_texture && texture(opacity_texture, fs_in.texcoord).r < 1.0)
		discard;

	// Diffuse color
	geometry_diffuse = vec4(0.0f);
	if (has_diffuse_texture)
		geometry_diffuse = vec4(0.3, 0.3, 0.3, 1.0) + 0.7*texture(diffuse_texture, vec2(0.5,0.25));

	// Specular color
	geometry_specular = vec4(0.0f);
	if (has_specular_texture)
		geometry_specular = vec4(0.0, 0.0, 0.0, 1.0);


	//texture normals range [0,1]
	geometry_normal = vec4(fs_in.normal * 0.5 + 0.5, 1.0);
	geometry_binormal = vec4(fs_in.binormal* 0.5 + 0.5, 1.0);
	geometry_tangent = vec4(fs_in.tangent * 0.5 + 0.5, 1.0);
}
