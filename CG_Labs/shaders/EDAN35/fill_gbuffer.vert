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

uniform mat4 vertex_model_to_world;

layout (location = 0) in vec3 vertex;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec3 texcoord;
layout (location = 3) in vec3 tangent;
layout (location = 4) in vec3 binormal;

out VS_OUT {
	vec4 world_pos;
	vec3 normal;
	vec2 texcoord;
	vec3 tangent;
	vec3 binormal;
	
	float kappa_max;
	float kappa_min;
	vec3 d_max;
	vec3 d_min;
} vs_out;


void main() {
	vs_out.world_pos = vertex_model_to_world * vec4(vertex, 1.0);
	
	vs_out.normal   = normalize(normal);
	vs_out.texcoord = texcoord.xy;
	vs_out.binormal = normalize(binormal);
	vs_out.tangent = normalize(tangent);
	//vs_out.tangent  = normalize(cross(normalize(binormal), normalize(normal)));

	vec3 normal = normalize(mat3(vertex_model_to_world) *normal);	
	vec3 tangent = normalize(mat3(vertex_model_to_world) *tangent);
	vec3 binormal = normalize(mat3(vertex_model_to_world) * binormal);

	//ambitions were high...
	mat3 curvature_tensor = mat3(1,0,0,0,1,0,0,0,1);
							//mat3(dot(tangent, tangent),	dot(tangent, binormal),		dot(tangent, normal),
							//	dot(binormal, tangent),	dot(binormal, binormal),	dot(binormal, normal),
							//	dot(normal, tangent),	dot(normal, binormal),		dot(normal, normal)
							//);
	
	float trace = 3;
				//curvature_tensor[0][0] + curvature_tensor[1][1] + curvature_tensor[2][2];
	float det = 1;
				//curvature_tensor[0][0] * (curvature_tensor[1][1] * curvature_tensor[2][2] - curvature_tensor[1][2] * curvature_tensor[2][1]) -
    			//curvature_tensor[0][1] * (curvature_tensor[1][0] * curvature_tensor[2][2] - curvature_tensor[1][2] * curvature_tensor[2][0]) +
    			//curvature_tensor[0][2] * (curvature_tensor[1][0] * curvature_tensor[2][1] - curvature_tensor[1][1] * curvature_tensor[2][0]);
	float disc = sqrt(max(0.0, trace * trace - 4.0 * det));
	float k_max = 0.5 * (trace + disc); //length(cross(normal,tangent));
	float k_min = 0.5 * (trace - disc); //length(cross(normal,binormal));

	//approx directions
	vec3 d_max = tangent;
	vec3 d_min = binormal;

	vs_out.kappa_max = k_max; // the size of the max "slope" (not really, were doing some serious simplifications here)
	vs_out.kappa_min = k_min;
	vs_out.d_max = d_max; //the dir 
	vs_out.d_min = d_min;
	

	gl_Position = camera.view_projection * vertex_model_to_world * vec4(vertex, 1.0);
}
