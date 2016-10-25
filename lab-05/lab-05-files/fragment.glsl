uniform sampler2D tex;
uniform sampler2D tex_am;
varying vec3 N;
varying vec4 position;
void main()
{
    vec4 color = texture2D(tex, gl_TexCoord[0].st);
    gl_FragColor = color;
    // vec4 color_am = texture2D(tex_am, gl_TexCoord[1].st);
    // gl_FragColor = mix(color, color_am, 1);
}
