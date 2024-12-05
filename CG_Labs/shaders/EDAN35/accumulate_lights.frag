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

layout (std140) uniform LightViewProjTransforms
{
	ViewProjTransforms lights[4];
};

uniform int light_index;

uniform sampler2D depth_texture;
uniform sampler2D normal_texture;
uniform sampler2D shadow_texture;

uniform vec2 inverse_screen_resolution;

uniform vec3 camera_position;

uniform vec3 light_color;
uniform vec3 light_position;
uniform vec3 light_direction;
uniform float light_intensity;
uniform float light_angle_falloff;

layout (location = 0) out vec4 light_diffuse_contribution;
layout (location = 1) out vec4 light_specular_contribution;

float toon(vec3 n, vec3 L) {
	float cos_a = max(0.0, dot(n,L));
	float toon_step = 3.0;
	float toon_shade = floor(cos_a * toon_step) / toon_step;
	return toon_shade;
}

void main()
{
	//texture space
	vec2 shadowmap_texel_size = 1.0f / textureSize(shadow_texture, 0);
	// screen space -> texture space [0,1]
	vec2 texcoords = vec2(gl_FragCoord.x, gl_FragCoord.y) * inverse_screen_resolution;
	
	//extract normal from texture, texture space. [0,1] -> [-1,1]
	vec3 normal = 2.0*vec3(texture(normal_texture, texcoords)) - vec3(1.0);
	normal = normalize(normal);
	//extract depth, texture space
	float depth = texture(depth_texture, texcoords).x;

	//2d -> 3d //screen space coord -> clip space coord, [0,1] -> [-1,1]
	vec4 projection_pos = vec4(texcoords, depth, 1.0) * 2 - 1; 
	// clip space -> homogeneous world space
	vec4 world_pos = camera.view_projection_inverse * projection_pos; 
	//1/w (perspective divide): homogeneous coord -> 3d coord
	//world space
	vec3 world_position = world_pos.xyz / world_pos.w;

	//shadows
	//world space -> shadow map space / light space
	vec4 light_space_pos = lights[light_index].view_projection * vec4(world_position,1.0);
	//perspective divide, homogeneous -> [-1,1]
	vec3 shadow_pixel = light_space_pos.xyz / light_space_pos.w;
	//[-1,1] -> [0,1], to align with texture coords
	//depth from light source
	float pixel_depth = shadow_pixel.z * 0.5 + 0.5;
	vec2 shadow_coords = shadow_pixel.xy * 0.5 + 0.5;
	//sample shadow map. [0,1]
	vec4 shadow_map_sample = texture(shadow_texture, shadow_coords);


	//smoother shadows
	int size = 5;
	//floor
	int range = size / 2;
	vec2 sample_coords;
	float sample_depth;
	// nbr of samples = 0%
	float shadow_intensity = size*size; 
	float epsilon = 0.00001;

	for(int i = -range; i <= range; ++i) {
		for(int j = -range; j <= range; ++j) {
			vec2 offset = vec2(i,j) * shadowmap_texel_size.xy;
			sample_coords = shadow_coords + offset;
			sample_depth = texture(shadow_texture,sample_coords).x;

			//depth comparison, if shadow
			if(sample_depth + epsilon < pixel_depth) {
				shadow_intensity -= 1.0;
			}
		}
	}
	//[0,1]
	float shadow_ratio = shadow_intensity / (size*size);

	//phong
	vec3 dist_vec = light_position - world_position.xyz; 
	vec3 L = normalize(dist_vec);

	vec3 V = normalize(camera_position - world_position.xyz);
	float diffuse = max(dot(normal , L), 0);
	float spec = pow(max(dot(reflect(-L, normal), V), 0), 1);


	//distance falloff
	float dist_sq = 1/dot(dist_vec, dist_vec);

	//angular falloff
	float theta = acos(dot(-L, normalize(light_direction)));
	float ang_falloff = max(0.0f, 1 - (theta / light_angle_falloff));

	float falloff = light_intensity * dist_sq * ang_falloff;

	float toon_color = toon(normal, L);

	vec3 col = toon_color * light_color;


	light_diffuse_contribution  = vec4(col, 1.0); //vec4(diffuse) * vec4(light_color, 1.0) * falloff * shadow_ratio ;
	light_specular_contribution = vec4(spec) * vec4(light_color, 1.0) * falloff * shadow_ratio;
}
