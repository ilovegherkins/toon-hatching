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
uniform sampler2D tangent_texture;
uniform sampler2D binormal_texture;
uniform sampler2D curvature_texture; // not working yet

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
uniform bool has_curve_hatching;

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

vec3 smooth_tangent(vec2 coords) {
	vec2 tex_size = 1.0/ textureSize(tangent_texture,0);
	vec3 smooth_tangent = vec3(0.0);
	for (int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			smooth_tangent += texture(tangent_texture, coords + vec2(i,j) *tex_size).xyz;
		}
	}
	smooth_tangent = normalize(smooth_tangent);
	return smooth_tangent;
}

float tangent_hatch(vec2 coords, vec3 world_pos, float color, int toon_step) {
	vec3 dir = smooth_tangent(coords);
	float proj = dot(world_pos, dir);

	//didnt need to rewrite this, we only use 3 steps anyway
	float hatch_color = 1.0;
	int shadow_level = int(ceil((1-color) * toon_step));
    for (int i = 0; i < toon_step; i++) {
        if (shadow_level > i) {
            float hatch_layer = step(fract(proj * hatch_spacing), 0.5);
            hatch_color *= mix(1.0, hatch_layer, float(i + 1) / float(toon_step)); 
        } else {
            break;
        }
    }
	return hatch_color;
}

//power iteration 
vec2 power(mat2 M, int it) {
    vec2 b_n = vec2(0.5, 0.4); // should be random
    for (int i = 0; i < it; ++i) {
        vec2 b_n1 = M * b_n;
        float b_n1_norm = length(b_n1);

        b_n = b_n1 / b_n1_norm;
    }
	return b_n;
}

vec3 get_curv_n_dir(vec2 coords){
	//partial derivatives of N with regards to t and b
    vec3 curr_n = texture(normal_texture, coords).xyz * 2.0 - 1.0; //should we sample or use the incoming? :think
    vec3 b = texture(binormal_texture, coords).xyz * 2.0 - 1.0;
    vec3 t = texture(tangent_texture, coords).xyz * 2.0 - 1.0;

    //sample normal along t & b plane
    float epsilon = 0.001;
    vec3 n_b = texture(normal_texture, coords + b.xy * epsilon).xyz * 2.0 - 1.0;
    vec3 n_t = texture(normal_texture, coords + t.xy * epsilon).xyz * 2.0 - 1.0;

    //approx derivatives
    vec3 d_ny = (n_b - curr_n)/epsilon;
    vec3 d_nx = (n_t - curr_n)/epsilon; 

    //shape operator
    mat2 S = mat2(dot(d_nx , t), dot(d_ny, t),
                dot(d_nx, b), dot(d_ny, b));


    //eigenvector approx
    vec2 p_dir = power(S, 5);

    //curvature is given by biggest eigenvalue
    float max_curv = max(p_dir.x, p_dir.y);
	
    return vec3(max_curv, p_dir);
}


float surface_hatch(vec2 coords, int steps) {
	vec3 curv_dir = get_curv_n_dir(coords);
    float max_curv = curv_dir.x;
    vec2 p_dir = normalize(curv_dir.yz);

	//restrict hatching 
	float min_d = 0.0;
	float max_d = 0.5;
	float d_scale = 0.3;
	//more stable than clamp
    float d = mix(min_d, max_d, smoothstep(0.0, 1.0, abs(max_curv) * d_scale));

    vec2 hatch_coords = coords.xy * p_dir;
    float hatch_lines = step(fract(hatch_coords.x * d), 0.5);

    float hatch = smoothstep(0.45, 0.55, hatch_lines); 

	return 1.0 - hatch;
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
	float curve_hatch = 1.0;

	if (has_hatching) {
		toon_color = toon(normal, L, falloff * shadow_ratio, toon_step);

		hatch_color = has_surf_hatching ? tangent_hatch(texcoords, world_position, toon_color, toon_step) : hatch(toon_color, toon_step);

		toon_color = 1.0;
	}
	if (has_toon) {
		toon_color = toon(normal, L, falloff * shadow_ratio, toon_step);
	}
	if (has_edges) {
		e = edge(texcoords, normal_texture, 1);
		e *= edge(texcoords, depth_texture, 1000);
	}
	if (has_curve_hatching) {
		curve_hatch = surface_hatch(texcoords, toon_step);
	}



	vec3 final_res = vec3(1.0) * toon_color * hatch_color * e * curve_hatch;

	light_diffuse_contribution  = vec4(vec3(final_res), 1.0); 
	light_specular_contribution = vec4(0.0); 
}
