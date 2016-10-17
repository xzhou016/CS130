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

    vec(const T& a, const T& b, const T& c,const T& d)
    {assert(n == 4);x[0]=a;x[1]=b;x[2]=c;x[3]= d;}

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
typedef vec<float,4> vec4;
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

