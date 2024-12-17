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
uniform sampler2D diffuse_texture;
uniform sampler2D direction_texture;

uniform vec2 inverse_screen_resolution;

uniform vec3 camera_position;

uniform vec3 light_color;
uniform vec3 light_position;
uniform vec3 light_direction;
uniform float light_intensity;
uniform float light_angle_falloff;

uniform float hatch_spacing;
uniform int hatch_sharpness;

layout (location = 0) out vec4 light_diffuse_contribution;
layout (location = 1) out vec4 light_specular_contribution;

float hatch(float color, int step) {
	int test = int(ceil(color * step));
	//hårdkodat, inte bra
	switch(test) {
		case 0:
			color = 0.0f;
			break;
		case 1:
			color *= abs(1-pow(sin(hatch_spacing*gl_FragCoord.x + hatch_spacing*gl_FragCoord.y), hatch_sharpness));
		case 2:
			color *= abs(1-pow(sin(hatch_spacing*gl_FragCoord.x - hatch_spacing*gl_FragCoord.y), hatch_sharpness));
		case 3:

			break;
	}
	return color;
}

vec3 hatch2(vec2 coords, vec3 world_pos, int toon_step, float toon_color) {
	vec3 surface_dir_sample = texture(direction_texture, coords).xyz * 2.0 - 1.0;
	surface_dir_sample = normalize(surface_dir_sample);

	float scale_factor = 0.9;
	float proj = dot(world_pos * scale_factor, surface_dir_sample);
	
	float hatch_pattern = abs(sin(proj * scale_factor * hatch_spacing * 3.14159)); 

	//test, doesnt look nice at all:)
	//float light_factor = clamp(1.0 - (spec), 0.0, 1.0);
	//hatch_pattern *= light_factor * hatch_mask;

	int hatch_level = int(ceil((1- toon_color) * toon_step));
    float hatch_mask = 0.0;

    if (hatch_level == 1) {
        hatch_mask = 1.0 - pow(hatch_pattern, 2); 
    } else if (hatch_level == 2) {
        hatch_mask = pow(hatch_pattern, 2); 
    }
	
	hatch_pattern = pow(hatch_pattern, float(hatch_sharpness));
	
	vec3 hatch_color = vec3(hatch_pattern);
	return hatch_color * hatch_mask;
}

float toon(vec3 n, vec3 L, float factor, int toon_step) {
	float cos_a = max(0.0, dot(n,L));
	float toon_shade = round(cos_a * toon_step * min(factor, 1.0)) / toon_step;
	return toon_shade;
}

void sobel_texture(inout vec4 n[9], sampler2D tex, vec2 coord) {
	float s_x = inverse_screen_resolution.x;
	float s_y = inverse_screen_resolution.y;
	n[0] = texture2D(tex, coord + vec2(-s_x, -s_y));
	n[1] = texture2D(tex, coord + vec2( 0.0, -s_y));
	n[2] = texture2D(tex, coord + vec2( s_x, -s_y));
	n[3] = texture2D(tex, coord + vec2(-s_x,  0.0));
	n[4] = texture2D(tex, coord);
	n[5] = texture2D(tex, coord + vec2( s_x,  0.0));
	n[6] = texture2D(tex, coord + vec2(-s_x,  s_y));
	n[7] = texture2D(tex, coord + vec2( 0.0,  s_y));
	n[8] = texture2D(tex, coord + vec2( s_x,  s_y));
}

float edge(vec2 t, sampler2D n_texture, sampler2D d_texture) {
	mat3 sobel_x = mat3( 1, 2, 1,
						 0, 0, 0,
						-1,-2,-1);
	mat3 sobel_y = transpose(sobel_x);

	vec2 offset;
	vec2 sample_pixel;
	float test_x = 0;
	float test_y = 0;
	for(int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			offset = vec2(i, j) * inverse_screen_resolution;
			sample_pixel = t + offset;
			test_x += length(texture(n_texture, sample_pixel).xyz) * sobel_x[i+1][j+1];
			test_x += length(texture(d_texture, sample_pixel).xyz) * sobel_x[i+1][j+1];
			test_y += length(texture(n_texture, sample_pixel).xyz) * sobel_y[i+1][j+1];
			test_y += length(texture(d_texture, sample_pixel).xyz) * sobel_y[i+1][j+1];
		}
	}
	float g = abs(test_x) + abs(test_y);
	g = step(0.5,g);

	return max(1-5*g, 0);
	
}

vec3 handle_curvature(vec2 coords, vec3 world_pos) {
	vec4 curv = texture(direction_texture, coords) * 2.0 - 1.0;
	float curv_w = curv.w;
	
	float proj = dot(world_pos, curv.xyz);
	float hatch = abs(1-pow(sin(proj * hatch_spacing), hatch_sharpness));
	hatch = hatch/(1.0 + abs(curv_w));

	float test_spacing = hatch_spacing * (1 + curv_w);
	float hatch_2 = abs(1-pow(sin(proj * test_spacing), hatch_sharpness));

	int scale = 300;
	//float hatch = abs(1-pow(sin(scale*hatch_spacing*curv.x + scale*hatch_spacing*curv.y + scale*hatch_spacing*curv.z), hatch_sharpness));

	float hatch_color = mix(hatch, hatch_2, 0.5);
	return hatch_color * vec3(1.0);

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
	int size = 7;
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
			sample_depth = texture(shadow_texture, sample_coords).x;

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

	const int toon_step = 3;

	float toon_color = toon(normal, L, falloff * shadow_ratio, toon_step);
	float hatch_color = hatch(toon_color, toon_step);
	vec3 hatch_color2 = hatch2(texcoords, world_position, toon_step, toon_color);
	hatch_color2 = vec3(toon_color) - hatch_color2;
	float e = edge(texcoords, normal_texture, depth_texture);
	vec3 final_res = vec3(e);

	vec3 in_from_vert = handle_curvature(texcoords, world_position);
	in_from_vert = in_from_vert* toon_color;

	vec3 debug = 1 - abs(world_position);

	light_diffuse_contribution  = vec4(in_from_vert, 1.0);//vec4(diffuse) * vec4(light_color, 1.0) * falloff * shadow_ratio;
	light_specular_contribution = vec4(0.0);//vec4(spec) * vec4(light_color, 1.0) * falloff * shadow_ratio;
}
