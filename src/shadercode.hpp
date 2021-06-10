
// the following is the shader code, which allows for the fuzzy look

namespace shader {
const char *lineVertex = R"(
#version 400

layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec4 vertexColor;
layout(location = 2) in vec2
    vertexSize;  // as a 2D texture cordinate, but we ignore the y

out Vertex {
  vec4 color;
  float size;
}
vertex;

uniform mat4 al_ModelViewMatrix;
uniform mat4 al_ProjectionMatrix;

void main() {
  gl_Position = al_ModelViewMatrix * vec4(vertexPosition, 1.0);
  vertex.color = vertexColor;
  vertex.size = vertexSize.x;
}
)";

const char *lineGeometry = R"(
#version 400

layout(lines) in;
layout(triangle_strip, max_vertices = 4) out;

uniform mat4 al_ProjectionMatrix;

in Vertex {
  vec4 color;
  float size;
}
vertex[];

out Fragment {
  vec4 color;
  float textureCoordinate;
}
fragment;

void main() {
  mat4 m = al_ProjectionMatrix;  // rename to make lines shorter
  vec4 a = gl_in[0].gl_Position;
  vec4 b = gl_in[1].gl_Position;

  const float r = 0.03;

  // does vec3(0.0, 0.0, 1.0) point at the eye in the coordinate system?
  // we hope that this billboards; it seems to, but really?
  // XXX i think that this is broken; it's just not quite right
  vec4 d = vec4(normalize(cross(b.xyz - a.xyz, vec3(0.0, 0.0, 1.0))), 0.0) * r;

  gl_Position = m * (a + d * vertex[0].size);
  fragment.color = vertex[0].color;
  fragment.textureCoordinate = 0.0;
  EmitVertex();

  gl_Position = m * (a - d * vertex[0].size);
  fragment.color = vertex[0].color;
  fragment.textureCoordinate = 1.0;
  EmitVertex();

  gl_Position = m * (b + d * vertex[1].size);
  fragment.color = vertex[1].color;
  fragment.textureCoordinate = 0.0;
  EmitVertex();

  gl_Position = m * (b - d * vertex[1].size);
  fragment.color = vertex[1].color;
  fragment.textureCoordinate = 1.0;
  EmitVertex();

  EndPrimitive();
}
)";

const char *lineFragment = R"(
#version 400

in Fragment {
  vec4 color;
  float textureCoordinate;
}
fragment;

uniform sampler1D alphaTexture;

layout(location = 0) out vec4 fragmentColor;

void main() {
  float a = texture(alphaTexture, fragment.textureCoordinate).r;
  if (a < 0.05) discard;
  fragmentColor = vec4(fragment.color.xyz, fragment.color.a * a);
}
)";

const char *pointVertex = R"(
#version 400

layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec4 vertexColor;
layout(location = 2) in vec2
    vertexSize;  // as a 2D texture cordinate, but we ignore the y

uniform mat4 al_ModelViewMatrix;
uniform mat4 al_ProjectionMatrix;

out Vertex {
  vec4 color;
  float size;
}
vertex;

void main() {
  gl_Position = al_ModelViewMatrix * vec4(vertexPosition, 1.0);
  vertex.color = vertexColor;
  vertex.size = vertexSize.x;
}
)";

const char *pointGeometry = R"(
#version 400

// take in a point and output a triangle strip with 4 vertices (aka a "quad")
//
layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

uniform mat4 al_ProjectionMatrix;

in Vertex {
  vec4 color;
  float size;
}
vertex[];

out Fragment {
  vec4 color;
  vec2 textureCoordinate;
}
fragment;

void main() {
  mat4 m = al_ProjectionMatrix;   // rename to make lines shorter
  vec4 v = gl_in[0].gl_Position;  // al_ModelViewMatrix * gl_Position

  // http://ogldev.atspace.co.uk/www/tutorial27/tutorial27.html
  //
  // in the tutorial at the link above, the billboarding computations
  // are done in worldspace using the camera position, but in this
  // example we apply the modelview matrix and do the billboarding
  // computations with respect to the origin which should be the
  // eye/camera position. are we wrong?

  float r = 0.15;
  r *= vertex[0].size;

  gl_Position = m * (v + vec4(-r, -r, 0.0, 0.0));
  fragment.textureCoordinate = vec2(0.0, 0.0);
  fragment.color = vertex[0].color;
  EmitVertex();

  gl_Position = m * (v + vec4(r, -r, 0.0, 0.0));
  fragment.textureCoordinate = vec2(1.0, 0.0);
  fragment.color = vertex[0].color;
  EmitVertex();

  gl_Position = m * (v + vec4(-r, r, 0.0, 0.0));
  fragment.textureCoordinate = vec2(0.0, 1.0);
  fragment.color = vertex[0].color;
  EmitVertex();

  gl_Position = m * (v + vec4(r, r, 0.0, 0.0));
  fragment.textureCoordinate = vec2(1.0, 1.0);
  fragment.color = vertex[0].color;
  EmitVertex();

  EndPrimitive();
}
)";

const char *pointFragment = R"(
#version 400

in Fragment {
  vec4 color;
  vec2 textureCoordinate;
}
fragment;

uniform sampler2D alphaTexture;

layout(location = 0) out vec4 fragmentColor;

void main() {
  float a = texture(alphaTexture, fragment.textureCoordinate).r;
  // if (a < 0.05) discard;
  fragmentColor = vec4(fragment.color.xyz, a);
}
)";
}  // namespace shader