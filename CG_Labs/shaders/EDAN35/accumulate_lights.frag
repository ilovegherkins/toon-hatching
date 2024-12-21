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
uniform bool has_hatching;
uniform bool has_toon;
uniform bool has_edges;
uniform bool has_surf_hatching;

layout (location = 0) out vec4 light_diffuse_contribution;
layout (location = 1) out vec4 light_specular_contribution;

float hatch(float color, int step) {
	int test = int(ceil(color * step));
	float res = has_toon ? color : 1.0;
	switch(test) {
		case 0:
		    res = 0.0f;
		    break;
		case 1:
		    res *= 1-pow(abs(sin(hatch_spacing*gl_FragCoord.x + hatch_spacing*gl_FragCoord.y)), hatch_sharpness);
		case 2:
		    res *= 1-pow(abs(sin(hatch_spacing*gl_FragCoord.x - hatch_spacing*gl_FragCoord.y)), hatch_sharpness);
		case 3:
		    break;
	}
	return res;
}

vec3 hatch2(vec2 coords, vec3 world_pos, int toon_step, float toon_color) {
	vec3 surface_dir_sample = texture(direction_texture, coords).xyz * 2.0 - 1.0;
	surface_dir_sample = normalize(surface_dir_sample);

	float scale_factor = 0.9;
	float proj = dot(world_pos * scale_factor, surface_dir_sample);
	
	float hatch_pattern = abs(sin(proj * scale_factor * hatch_spacing * 3.14159)); 

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

float edge(vec2 t, sampler2D tex, int weight) {
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
			test_x += length(texture(tex, sample_pixel).xyz) * sobel_x[i+1][j+1] * weight;
			test_y += length(texture(tex, sample_pixel).xyz) * sobel_y[i+1][j+1] * weight;
		}
	}
	float g = abs(test_x) + abs(test_y);
	g = step(0.5,g);

	return max(1-5*g, 0);
	
}

vec3 smooth_dir(vec2 coords) {
	vec2 tex_size = 1.0/ textureSize(direction_texture,0);
	vec3 smooth_dir = vec3(0.0);
	for (int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			smooth_dir += texture(direction_texture, coords + vec2(i,j) *tex_size).xyz;
		}
	}
	smooth_dir = normalize(smooth_dir);
	return smooth_dir;
}

vec3 rotate_dir(vec3 dir, float ang) {
	float angle = radians(ang);
	float c = cos(angle);
	float s = sin(angle);
	return vec3(
		c * dir.x - s * dir.y,
		s * dir.x + c * dir.y,
		dir.z
	);
}

float surface_hatch(vec2 coords, vec3 world_pos, int method, float color, int step) {
	vec3 dir = vec3(0.0);
	vec4 surf = texture(direction_texture, coords) * 2.0 - 1.0; 
	

	switch (method) {
		case 0:
			// flowerpots look ok here and curtains look below decent
			dir = smooth_dir(coords); 
			break;
		case 1:
			// flowerpots look good here, but curtains are an extreme mess
			dir = surf.xyz;
			break;
		default:
			//defaulting to the one that looks less bad
			dir = smooth_dir(coords); 
			break;
	}
	float k_max = surf.w; // max surf 

	//flat surfaces leads to wide spacing of hatching lines (does nothing rn, k_max = 1.0)
	float surface_bend = abs(k_max);
	float dyn_spacing = hatch_spacing * surface_bend;
	//dir = length(dir) > 0.1 ? dir : vec3(0.0, 1.0, 0.0); 

	float proj1 = dot(world_pos, dir);
	float proj2 = dot(world_pos, rotate_dir(dir, 45.0));

	float hatch_color = 1.0;
	
	int shadow_level = int(ceil((1-color) * step));

    for (int i = 0; i < step; i++) {
        if (shadow_level > i) {
            //float proj = (i % 2 == 0) ? proj1 : proj2; // alternated dir depending on shadow depth, doesnt look good with current impl.
			float proj = proj1; 
            float hatch_layer = 1.0 - pow(abs(sin(proj * dyn_spacing)), hatch_sharpness);
            hatch_color *= mix(1.0, hatch_layer, float(i + 1) / float(step)); 
        } else {
            break;
        }
    }

	return hatch_color;
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
			vec2 offset = vec2(i,j) * shadowmap_texel_size;
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
	float toon_color = 1.0;
	float hatch_color = 1.0;
	float e = 1.0;

	int surf_hatch_method = 0; // 0 or 1. 0 looks less bad in general
	if (has_hatching) {
		toon_color = toon(normal, L, falloff * shadow_ratio, toon_step);

		hatch_color = has_surf_hatching ? surface_hatch(texcoords, world_position, surf_hatch_method, toon_color, toon_step) : hatch(toon_color, toon_step);

		toon_color = 1.0;
	}
	if (has_toon) {
		toon_color = toon(normal, L, falloff * shadow_ratio, toon_step);
	}
	if (has_edges) {
		e = edge(texcoords, normal_texture, 1);
		e *= edge(texcoords, depth_texture, 1000);
	}

	vec3 final_res = vec3(1.0) * toon_color * hatch_color * e;
	light_diffuse_contribution  = vec4(final_res, 1.0); //vec4(diffuse) * vec4(light_color, 1.0) * falloff * shadow_ratio;
	light_specular_contribution = vec4(0.0); //vec4(spec) * vec4(light_color, 1.0) * falloff * shadow_ratio;
}
