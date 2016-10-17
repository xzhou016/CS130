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
 vector<mglVertices> vertex_set;

struct mglMatrix {
     MGLfloat fMatrix[16];

     mglMatrix(){make_zero();}

     void make_zero()
    {for(size_t i = 0; i < 16; i++) fMatrix[i] = 0;}

    void make_id()
    {make_zero();for(int i = 0; i < 16; i+=5) fMatrix[i] = 1;}

    mglMatrix& operator = (const MGLfloat *matrix){
        int other_mat_count = 0;
        for (size_t i = 0; i < 16; i++) {
                fMatrix[i] = matrix[other_mat_count];
                other_mat_count++;
        }
        return *this;
    }

};
mglMatrix render_matrix;

stack<mglMatrix> render_stack;

 MGLbool render_started;
 MGLint  render_poly_mode;
 MGLint  render_matrix_state;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
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
    // assert(width > 0 && height > 0 ); //make sure w, h are positive
    //
    // //create 2-D display
    // MGLpixel screen[width*height];
    // MGLfloat z_buffer[width*height];
    //
    // for (size_t i = 0; i < width*height; i++) {
    //     screen[i]   = data[i]; //read raw data, should be 0
    //     z_buffer[i] = -1;      //z is at infinity
    // }
    //
    // render_matrix.make_id();


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
    // switch (render_poly_mode) {
    //     case MGL_TRIANGLES: {
    //
    //         break;
    //     }
    //     case MGL_QUADS: {
    //         break;
    //     }
    //
    //     default :
    //         MGL_ERROR("Incompelete render!");
    // }
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
    mglVertices node;
    node = {.vert_x = x, .vert_y = y, .vert_z = x, .vert_w = 1};

    vertex_set.push_back(node);

}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    render_matrix_state = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
    // render_stack.push(render_matrix);
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    // render_matrix = render_stack.top();
    // render_stack.pop();
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
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
void mglLoadMatrix(const MGLfloat *matrix)
{
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
void mglMultMatrix(const MGLfloat *matrix)
{
    mglMatrix   mat_sum;
    MGLsize     sum_index = 0;

    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            mat_sum.fMatrix[sum_index] = matrix[j] * render_matrix.fMatrix[i]
                                + matrix[j+4] * render_matrix.fMatrix[i + 1]
                                + matrix[j+8] * render_matrix.fMatrix[i + 2]
                                + matrix[j+12] * render_matrix.fMatrix[i + 3];
            sum_index++;
        }
    }

    render_matrix = mat_sum;
}

// void mglMultMatrix(MGLfloat x,
//               MGLfloat y,
//               MGLfloat z)
// {
//     mglMatrix   mat_sum;
//     MGLsize     sum_index = 0;
//
//     for (size_t i = 0; i < 4; i++) {
//         for (size_t j = 0; j < 4; j++) {
//             mat_sum.fMatrix[sum_index] = matrix[j] * render_matrix.fMatrix[i]
//                                 + matrix[j+4] * render_matrix.fMatrix[i + 1]
//                                 + matrix[j+8] * render_matrix.fMatrix[i + 2]
//                                 + matrix[j+12] * render_matrix.fMatrix[i + 3];
//             sum_index++;
//         }
//     }
//
//     render_matrix = mat_sum;
// }

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
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
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{

}
