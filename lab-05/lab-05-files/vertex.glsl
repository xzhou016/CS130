uniform sampler2D tex;
varying vec3 N;
varying vec4 position;

void main()
{
    gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0;
    gl_Position = ftransform();
}
