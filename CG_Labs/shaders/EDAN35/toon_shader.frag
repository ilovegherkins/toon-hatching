#version 410

uniform vec3 light_position;
uniform vec3 light_direction;

out vec4 FragColor;

int main() {
    FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}