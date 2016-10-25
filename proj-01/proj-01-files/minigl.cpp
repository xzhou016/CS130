/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 */

#include "minigl.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <iostream>
#include <stack>


using namespace std;

/**
 * Global values for rasterization
 */
struct mglVertices {
    MGLfloat vert_x;
    MGLfloat vert_y;
    MGLfloat vert_z;
    MGLfloat vert_w;

    MGLpixel pixel_color;
};

//global vector containing all processed vertices
vector<mglVertices> vertex_set;

struct mglMatrix {
     MGLfloat fMatrix[16];

     mglMatrix(){make_zero();}

     //= the zero matrix
     void make_zero(){
         for(size_t i = 0; i < 16; i++)
            fMatrix[i] = 0;
     }
    //= the identity matrix
    void make_id(){
        make_zero();
        for(int i = 0; i < 16; i+=5)
            fMatrix[i] = 1;
    }

    mglMatrix& operator = (const MGLfloat *matrix){
        int other_mat_count = 0;
        for (size_t i = 0; i < 16; i++) {
                fMatrix[i] = matrix[other_mat_count];
                other_mat_count++;
        }
        return *this;
    }

    //for multiplying local matrix with a given matrix
    mglMatrix& operator * (mglMatrix B){
        MGLfloat sum[16];

        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
                sum[j*4 + i] = B.fMatrix[i]     *fMatrix[j*4]
                            + B.fMatrix[i + 4]  *fMatrix[j*4 + 1]
                            + B.fMatrix[i + 8]  *fMatrix[j*4 + 2]
                            + B.fMatrix[i + 12] *fMatrix[j*4 + 3];
            }
        }

        for (size_t i = 0; i < 16; i++) {
            fMatrix[i] = sum[i];
        }
        return *this;
    }

    //for multiplying vertices with the local matrix
    mglVertices multVerex(mglVertices in_vertex){
        mglVertices node;
        vector<MGLfloat> v;
        for (size_t i = 0; i < 4; i++) {
            v.push_back(  fMatrix[i]    *in_vertex.vert_x
                        + fMatrix[i+4]  *in_vertex.vert_y
                        + fMatrix[i+8]  *in_vertex.vert_z
                        + fMatrix[i+12] *in_vertex.vert_w);
        }
        node = {  .vert_x = v.at(0)
                , .vert_y = v.at(1)
                , .vert_z = v.at(2)
                , .vert_w = v.at(3)};
        return node;
    }

};

/*******************
 *Global variables
 */
mglMatrix render_matrix;   //4x4 matrix for transformation

stack<mglMatrix> render_model_stack;   //model view
stack<mglMatrix> render_proj_stack;    //projection

MGLbool        render_started;               //check if render started
MGLint         render_poly_mode;             //choose between v3 and quad
MGLmatrix_mode render_state = MGL_MODELVIEW; //FSM for choosing modelview and projection
MGLpixel*      render_bary_color;
MGLpixel*      render_z_buff;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


MGLfloat comput_area(mglVertices a, mglVertices b, mglVertices c){
    return (b.vert_x - a.vert_x)*(c.vert_y - a.vert_y) - (b.vert_y - a.vert_y)*(c.vert_x - a.vert_x);
}
/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    MGLpixel color = 000000;

    mglVertices a;
    mglVertices b;
    mglVertices c;
    mglVertices p;
    MGLfloat tri_area = comput_area(a,b,c);

    for (size_t i = 0; i < width; i++) {
        p.vert_x = i;
        for (size_t j = 0; j < height; j++) {
            p.vert_y = j;

            MGLfloat alpha = comput_area(p,b,c)/tri_area;
            MGLfloat beta = comput_area(a,p,c)/tri_area;
            MGLfloat gamma = comput_area(a,b,p)/tri_area;

            MGLfloat interpol_c = alpha*a.vert_z + beta*b.vert_z + gamma*c.vert_z;
            data[width*j + i] = color
        }
    }
    //push it into the framebuffer

}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    render_started = true;
    render_poly_mode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{

}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    mglVertex3(x, y, 0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    mglVertices node; //create a vertex
    node = {.vert_x = x, .vert_y = y, .vert_z = x, .vert_w = 1};

    //get the modelview * projection matrix
    mglMatrix transformed_geometry =  render_model_stack.top()
                                    * render_proj_stack.top();

    //do proj*modelview*vertex
    vertex_set.push_back(transformed_geometry.multVerex(node));
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode){
    switch (mode) {
        case MGL_MODELVIEW: {
            if(render_state != mode) {
                render_state = MGL_PROJECTION;
                mglPushMatrix();
            }
            //render_matrix = render_model_stack.top();
            render_state = mode;
            break;
        }
        case MGL_PROJECTION: {
            if (render_state != mode) {
                render_state = MGL_MODELVIEW;
                mglPushMatrix();
            }
            //render_matrix = render_proj_stack.top();
            render_state = mode;
            break;
        }
        default:    MGL_ERROR("Cannot find mode"); break;
    }

}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix(){
    switch (render_state) {
        case MGL_MODELVIEW:
            render_model_stack.push(render_matrix);break;
        case MGL_PROJECTION:
            render_proj_stack.push(render_matrix); break;
        default:
            MGL_ERROR("Cannot push onto stack");
    }
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix(){
    switch (render_state) {
        case MGL_MODELVIEW:{
            render_matrix = render_model_stack.top();
            render_model_stack.pop();
            break;
        }
        case MGL_PROJECTION: {
            render_matrix = render_proj_stack.top();
            render_proj_stack.pop();
            break;
        }

        default:
            MGL_ERROR("Cannot pop stack");
    }
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity(){
    render_matrix.make_id();
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix){
    render_matrix.make_zero();
    render_matrix = matrix;
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix){
    mglMatrix   mat_sum;
    MGLsize     sum_index = 0;

    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            mat_sum.fMatrix[sum_index] = matrix[j*4] * render_matrix.fMatrix[i]
                                        + matrix[j*4+4] * render_matrix.fMatrix[i + 4]
                                        + matrix[j*4+8] * render_matrix.fMatrix[i + 8]
                                        + matrix[j*4+12] * render_matrix.fMatrix[i + 12];
            sum_index++;
        }
    }

    render_matrix = mat_sum;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mglMatrix trans_matrix;
    trans_matrix.fMatrix[12] = x;
    trans_matrix.fMatrix[13] = y;
    trans_matrix.fMatrix[14] = z;

    mglMultMatrix(trans_matrix.fMatrix);

}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
    mglMatrix temp_rotate_matrix;
    MGLfloat rotate_cos = cos(angle*3.1415926/180);
    MGLfloat rotate_sin = sin(angle*3.1415926/180);

    //normalize
    MGLfloat magnitude = sqrt(x*x + y*y + z*z);
    x = x/magnitude;
    y = y/magnitude;
    z = z/magnitude;

    //make rotation from given axis and angle
    //column 1
    temp_rotate_matrix.fMatrix[0] = rotate_cos + x*x*(1 - rotate_cos);
    temp_rotate_matrix.fMatrix[1] = y*x*(1 - rotate_cos) + z*rotate_sin;
    temp_rotate_matrix.fMatrix[2] = z*x*(1 - rotate_cos) - y*rotate_sin;
    temp_rotate_matrix.fMatrix[3] = 0;
    //column 2
    temp_rotate_matrix.fMatrix[4] = x*x*(1 - rotate_cos) - z*rotate_sin;
    temp_rotate_matrix.fMatrix[5] = rotate_cos + y*y*(1-rotate_cos);
    temp_rotate_matrix.fMatrix[6] = z*y*(1 - rotate_cos) + x*rotate_sin;
    temp_rotate_matrix.fMatrix[7] = 0;
    //column 3
    temp_rotate_matrix.fMatrix[8] = x*z*(1 - rotate_cos) + y*rotate_sin;
    temp_rotate_matrix.fMatrix[9] = y*z*(1 - rotate_cos) -   x*rotate_sin;
    temp_rotate_matrix.fMatrix[10] = rotate_cos + z*z*(1-rotate_cos);
    temp_rotate_matrix.fMatrix[11] = 0;
    //column 4
    temp_rotate_matrix.fMatrix[12] = 0;
    temp_rotate_matrix.fMatrix[13] = 0;
    temp_rotate_matrix.fMatrix[14] = 0;
    temp_rotate_matrix.fMatrix[15] = 1;

    render_matrix = render_matrix*temp_rotate_matrix;
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mglMatrix scale_mat;
    scale_mat.fMatrix[0] = x;
    scale_mat.fMatrix[5] = y;
    scale_mat.fMatrix[10] = z;

    mglMultMatrix(scale_mat.fMatrix);
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
    mglMatrix view;
    view.make_id();

    view.fMatrix[0] = (2*near)/(right - left);
    view.fMatrix[5] = (2*near)/(top - bottom);

    view.fMatrix[8] = (right + left)/(right - left);
    view.fMatrix[9] = (top + bottom)/(top - bottom);
    view.fMatrix[10] = -1*(far + near)/(far - near);
    view.fMatrix[11] = -1;

    view.fMatrix[14] = -2*far*near/(far - near);
    view.fMatrix[15] = 0;

    render_matrix = render_matrix*view;
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    mglMatrix compute_ortho;

    compute_ortho.fMatrix[ 0] = 2 /(right - left);
    compute_ortho.fMatrix[ 5] = 2 /(top - bottom);
    compute_ortho.fMatrix[10] = -2/(far - near);
    compute_ortho.fMatrix[12] = -(right+left)/(right - left);
    compute_ortho.fMatrix[13] = -(top+bottom)/(top - bottom);
    compute_ortho.fMatrix[14] = -(far+near)  /(far - near);

    render_matrix = render_matrix*compute_ortho;
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{

}
