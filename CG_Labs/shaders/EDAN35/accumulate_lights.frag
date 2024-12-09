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
			test_y += length(texture(n_texture, sample_pixel).xyz) * sobel_y[i+1][j+1];
		}
	}
	float g = abs(test_x) + abs(test_y);

	return max(1-5*g, 0);

	//vec4[9] depth;
	//vec4[9] normal;
	//sobel_texture(depth, d_texture, t);
	//sobel_texture(normal, n_texture, t);
//
	//vec4 sobel_edge_depth_h = 1.0*depth[2] + (2.0*depth[5]) + 1.0*depth[8] - (1.0*depth[0] + (2.0*depth[3]) + 1.0*depth[6]);
  	//vec4 sobel_edge_depth_v = 1.0*depth[0] + (2.0*depth[1]) + 1.0*depth[2] - (1.0*depth[6] + (2.0*depth[7]) + 1.0*depth[8]);
//
	//vec4 sobel_edge_normal_h = 1.0*normal[2] + (2.0*normal[5]) + 1.0*normal[8] - (1.0*normal[0] + (2.0*normal[3]) + 1.0*normal[6]);
  	//vec4 sobel_edge_normal_v = 1.0*normal[0] + (2.0*normal[1]) + 1.0*normal[2] - (1.0*normal[6] + (2.0*normal[7]) + 1.0*normal[8]);
//
	//vec4 sobel = sqrt((sobel_edge_depth_h * sobel_edge_depth_h) + (sobel_edge_depth_v * sobel_edge_depth_v))
	//           ;//+ sqrt((sobel_edge_normal_h * sobel_edge_normal_h) + (sobel_edge_normal_v * sobel_edge_normal_v));
//
	//return 1.0-sobel.rgb;
	
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

	float e = edge(texcoords, normal_texture, depth_texture);
	
	vec3 col = e * hatch_color * vec3(1.0);
	
	light_diffuse_contribution  = vec4(col, 1.0);//vec4(diffuse) * vec4(light_color, 1.0) * falloff * shadow_ratio;
	light_specular_contribution = vec4(spec) * vec4(light_color, 1.0) * falloff * shadow_ratio;
}
