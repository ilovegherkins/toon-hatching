#version 410

uniform sampler2D normals_texture;
uniform sampler2D tangents_texture;
uniform sampler2D binormals_texture;

uniform mat4 normal_model_to_world; 
uniform mat4 vertex_model_to_world;

uniform vec2 inverse_screen_resolution;

in VS_OUT {
    //vec4 world_pos;
	vec3 normal;
	vec2 texcoord;
	vec3 tangent;
	vec3 binormal;

    //vec3 d_nx;
    //vec3 d_ny;
} fs_in;

layout (location = 0) out vec4 curvature;

//vec2 power(mat2 M, int it) {
//    vec2 b_n = vec2(1,4); // should be random
//    float b_n1 = 0.0;
//    float b_n1_norm = 0.0;
//    for (int i = 0; i < it; ++i) {
//        float b_n1 = dot(M, b_n);
//        b_n1_norm = length(b_n1);
//
//        b_n = b_n1 / b_n1_norm;
//    }
//}

void main() {
    //this isnt working rn and i dont have more time
    //vec2 coords = vec2(gl_FragCoord.x, gl_FragCoord.y) * inverse_screen_resolution;

    //partial derivatives of n with regards to t and b 
    //vec3 curr_n = texture(normals_texture, coords).xyz; //should we sample or use the incoming? :think
    //vec3 b = texture(binormals_texture, coords).xyz;
    //vec3 t = texture(tangents_texture, coords).xyz;

    //sample normal along t & b plane
    //float epsilon = 0.001;
    //vec3 n_b = texture(normal_texture, b.xy * epsilon).xyz;
    //vec3 n_t = texture(normal_texture, t.xy * epsilon).xyz;

    //approx derivatives
    //vec3 d_ny = (n_b - curr)/epsilon;
    //vec3 d_nx = (n_t - curr)/epsilon; 

    //shape operator
    //mat2 S = mat2(dot(d_nx ,fs_in.tangent), dot(d_ny, fs_in.tangent),
    //            dot(d_nx, fs_in.bitangent), dot(d_ny, fs_in.bitangent));

    //eigenvalue approx - power iteration 
    //vec2 p_dir = power(S, 5);

    //curvature is given by biggest eigenvalue
    //float max_curv = max(p_dir);

    //curvature = vec4(max_curv, p_dir, 1.0);
    //vec3 n = texture(normal_texture, texcoord);
    curvature = vec4(fs_in.normal,0.0);

}