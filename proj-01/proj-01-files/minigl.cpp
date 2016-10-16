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
 * vec.h
 **********************************/
template<class T, int n> struct vec;
template<class T, int n> T dot(const vec<T,n>& u,const vec<T,n>& v);

template<class T, int n>
struct vec
{
    T x[n];

    vec()
    {make_zero();}

    vec(const T& a)
    {assert(n == 1);x[0]=a;}

    vec(const T& a, const T& b)
    {assert(n == 2);x[0]=a;x[1]=b;}

    vec(const T& a, const T& b, const T& c)
    {assert(n == 3);x[0]=a;x[1]=b;x[2]=c;}

    void make_zero()
    {for(int i = 0; i < n; i++) x[i] = 0;}

    vec& operator += (const vec& v)
    {for(int i = 0; i < n; i++) x[i] += v.x[i]; return *this;}

    vec& operator -= (const vec& v)
    {for(int i = 0; i < n; i++) x[i] -= v.x[i]; return *this;}

    vec& operator *= (const vec& v)
    {for(int i = 0; i < n; i++) x[i] *= v.x[i]; return *this;}

    vec& operator /= (const vec& v)
    {for(int i = 0; i < n; i++) x[i] /= v.x[i]; return *this;}

    vec& operator *= (const T& c)
    {for(int i = 0; i < n; i++) x[i] *= c; return *this;}

    vec& operator /= (const T& c)
    {for(int i = 0; i < n; i++) x[i] /= c; return *this;}

    vec operator + () const
    {return *this;}

    vec operator - () const
    {vec r; for(int i = 0; i < n; i++) r[i] = -x[i]; return r;}

    vec operator + (const vec& v) const
    {vec r; for(int i = 0; i < n; i++) r[i] = x[i] + v.x[i]; return r;}

    vec operator - (const vec& v) const
    {vec r; for(int i = 0; i < n; i++) r[i] = x[i] - v.x[i]; return r;}

    vec operator * (const vec& v) const
    {vec r; for(int i = 0; i < n; i++) r[i] = x[i] * v.x[i]; return r;}

    vec operator / (const vec& v) const
    {vec r; for(int i = 0; i < n; i++) r[i] = x[i] / v.x[i]; return r;}

    vec operator * (const T& c) const
    {vec r; for(int i = 0; i < n; i++) r[i] = x[i] * c; return r;}

    vec operator / (const T& c) const
    {vec r; for(int i = 0; i < n; i++) r[i] = x[i] / c; return r;}

    const T& operator[] (int i) const
    {return x[i];}

    T& operator[] (int i)
    {return x[i];}

    T magnitude_squared() const
    {return dot(*this, *this);}

    T magnitude() const
    {return sqrt(magnitude_squared());}

    // Be careful to handle the zero vector gracefully
    vec normalized() const
    {T mag = magnitude(); if(mag) return *this / mag; vec r; r[0] = 1; return r;};
};

template <class T, int n>
vec<T,n> operator * (const T& c, const vec<T,n>& v)
{return v*c;}

template <class T, int n>
T dot(const vec<T,n> & u, const vec<T,n> & v)
{
    T r  =  0;
    for(int i = 0; i < n; i++) r += u.x[i] * v.x[i];
    return r;
}

template <class T >
vec<T,3> cross(const vec<T,3> & u, const vec<T,3> & v)
{
    return vec<T,3> (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0]);
}

template <class T, int n>
std::ostream& operator << (std::ostream& out, const vec<T,n> & u)
{
    for(int i = 0; i < n; i++)
    {
        if(i) out << ' ';
        out << u[i];
    }
    return out;
}

typedef vec<float,2> vec2;
typedef vec<float,3> vec3;
typedef vec<int,2> ivec2;



/**
 * mat.h
**********************************************/
template<class T, int m, int n> struct mat;
template<class T, int m, int n> T dot(const mat<T,m,n>& u,const mat<T,m,n>& v);

template<class T, int m, int n>
struct mat
{
    T x[m][n];

    mat()
    {make_zero();}

    void make_zero()
    {for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) x[i][j] = 0;}

    void make_id()
    {assert(m==n);make_zero();for(int i = 0; i < m; i++) x[i][i] = 1;}

    mat& operator = (const MGLfloat *matrix){
        int other_mat_count = 0;
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                x[i][j] = matrix[other_mat_count];
                other_mat_count++;
            }
        }
        return *this;
    }

    mat& operator += (const mat& v)
    {for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) x[i][j] += v.x[i][j]; return *this;}

    mat& operator -= (const mat& v)
    {for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) x[i][j] -= v.x[i][j]; return *this;}

    mat& operator *= (const T& c)
    {for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) x[i][j] *= c; return *this;}

    mat& operator /= (const T& c)
    {for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) x[i][j] /= c; return *this;}

    mat operator + () const
    {return *this;}

    mat operator - () const
    {mat r; for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) r[i] = -x[i][j]; return r;}

    mat operator + (const mat& v) const
    {mat r; for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) r[i] = x[i][j] + v.x[i][j]; return r;}

    mat operator - (const mat& v) const
    {mat r; for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) r[i] = x[i][j] - v.x[i][j]; return r;}

    mat operator * (const T& c) const
    {mat r; for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) r[i] = x[i][j] * c; return r;}

    mat operator / (const T& c) const
    {mat r; for(int i = 0; i < m; i++) for(int j = 0; j < n; j++) r[i] = x[i][j] / c; return r;}

    template<int p>
    mat<T,m,p> operator * (const mat<T,n,p>& M) const
    {
        mat<T,m,p> r;
        for(int i = 0; i < m; i++)
            for(int j = 0; j < n; j++)
                for(int k = 0; k < p; k++)
                    r.x[i][k] += x[i][j] * M.x[j][k];
        return r;
    }

    const T& operator() (int i, int j) const
    {return x[i][j];}

    T& operator() (int i, int j)
    {return x[i][j];}
};

template <class T, int m, int n>
mat<T,m,n> operator * (const T& c, const mat<T,m,n>& v)
{return v*c;}

template <class T, int m, int n>
std::ostream& operator << (std::ostream& out, const mat<T,m,n> & M)
{
    for(int i = 0; i < m; i++)
    {
        if(i) out << " ; ";
        for(int j = 0; j < n; j++)
        {
            if(j) out << ' ';
            out << M.x[i][j];
        }
    }
    return out;
}

typedef mat<float,4,4> mat4;

// 4x4 rotation matrix taking vector "from" into vector "to"
//
inline mat4 from_rotated_vector(const vec3& from,const vec3& to)
{
    vec3 A[3]={from.normalized()};
    vec3 B[3]={to.normalized()};
    A[1]=B[1]=cross(A[0],B[0]).normalized();
    A[2]=cross(A[0],A[1]);
    B[2]=cross(B[0],B[1]);

    mat4 M;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                M(i,j)+=B[k][i]*A[k][j];
    M(3,3)=1;
    return M;
}

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

/**
 * Global values for rasterization
 */
 MGLbool        render_started;
 mat4           render_matrix;
 vec3           render_vertex_3;
 vec2           render_vertex_2;
 MGLint         render_poly_state, render_matrix_state;
 stack<mat4>    render_stack;


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
    assert(width > 0 && height > 0 ); //make sure w, h are positive

    //create 2-D display
    MGLpixel screen[width*height];
    MGLfloat z_buffer[width*height];

    for (size_t i = 0; i < width*height; i++) {
        screen[i]   = data[i]; //read raw data, should be 0
        z_buffer[i] = -1;      //z is at infinity
    }

    mat4 raster_matrix;
    raster_matrix.make_id();


}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    render_started = true;
    render_poly_state = mode;
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

}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
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
    render_stack.push(render_matrix);
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    render_matrix = render_stack.top();
    render_stack.pop();
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

}

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
