/*! \file smtensor.h
\brief Templated classes to work with (fixed size) tensors.
*/
#ifndef SMTENSOR_H
#define SMTENSOR_H

#include <array>
#include <cassert>
#include <cmath>
#include "config.h"

namespace EOS_Toolkit {

enum zero_literal {ZERO=0}; 
enum one_literal {ONE=1};

template<class T, int N> class sm_matrix_sym;
template<class T, int N> class sm_matrix_sqr;
template<class T, int N, bool UP1, bool UP2> class sm_tensor2_sym;

namespace detail {
template<int N, class... A> 
using fix_num_args = typename std::enable_if<(sizeof...(A) == N) && 
                                             (sizeof...(A)>1)>::type;
}

//-------------------------------------------------------------------
//  Raw vector class
//-------------------------------------------------------------------

/**\brief Class representing fixed size vector

\tparam T Data type of components
\tparam N Dimensionality
**/
template<class T, int N> 
class sm_vector {
  static_assert(N>0, "Vector dimension must be > 0");
  
  using me_t = sm_vector<T,N>;
  
  std::array<T,N> v;
  
  public:
  
  using elem_t = T;
  
  enum {SIZE=N};
  
  sm_vector()                       noexcept = default;
  sm_vector(me_t&& that)            noexcept = default;
  sm_vector(const me_t& that)       noexcept = default;
  sm_vector& operator=(const me_t& that) noexcept = default;
  sm_vector& operator=(me_t&& that)      noexcept = default;
  

  template<class T2>
  explicit sm_vector(const sm_vector<T2,N>& that) noexcept 
  {
    for (int i=0; i<N; ++i) {
      v[i] = that(i); 
    }  
  }
  template<class T2>
  sm_vector& operator=(const sm_vector<T2,N>& that) noexcept 
  {
    *this = sm_vector<T,N>(that);
    return *this;
  }
  
  template<class... An, 
           class Na=detail::fix_num_args<N,An...> > 
  explicit sm_vector(An&&... an) 
  : v{{std::forward<An>(an)...}} {}
    
  sm_vector(zero_literal) {zero();}

  T& operator()(int j) {return v.at(j);}
  const T& operator()(int j) const {return v.at(j);}

  template<class T2>
  void operator*=(T2 z) {
    for(int i=0; i<N; ++i) {
      v[i] *= z;
    }
  }
  
  template<class T2>
  void operator/=(T2 z) {
    for(int i=0; i<N; ++i) {
      v[i] /= z;
    }
  }
  
  template<class T2>
  void operator+=(const sm_vector<T2, N>& a) {
    for(int i=0; i<N; ++i) {
      v[i] += a(i);
    }
  }
  
  template<class T2>
  void operator-=(const sm_vector<T2,N>& a) {
    for(int i=0; i<N; ++i) {
      v[i] -= a(i);
    }
  }

  void negate() {
    for(int i=0; i<N; ++i) {
      v[i] = -v[i];
    }
  }

  T dot(const me_t& a) const {
    T erg=v[0]*a(0); 
    for (int i=1; i<N; ++i) {
      erg += v[i]*a(i);
    }
    return erg;
  }

  void assign_prod(const sm_matrix_sym<T, N>& m, const me_t& w);
  void assign_prod(const me_t& v, const sm_matrix_sym<T, N>& m) {
    assign_prod(m,v);
  }
  void assign_prod(const sm_matrix_sqr<T, N>& m, const me_t& w);
  void assign_prod(const me_t& w, const sm_matrix_sqr<T, N>& m);


  me_t operator+(const me_t& a) const {
    me_t e{*this};
    e += a; 
    return e;
  }
  me_t operator-(const me_t& a) const {
    me_t e{*this};
    e -= a; 
    return e;
  }
  me_t operator*(T z) const {
    me_t e{*this}; 
    e *= z;
    return e;
  }
  me_t operator/(T z) const {
    me_t e{*this};
    e /= z; 
    return e;
  }
  me_t operator-() const {
    me_t e{*this};
    e.negate();
    return e;
  }

  void zero() {
    v.fill(0);
  }
};

template<class T, int N>
auto operator*(const T& a, const sm_vector<T, N> &v)
-> sm_vector<T, N>
{
  return v * a;
}


//--------------------------------------------------------------------
//  Raw symmetric matrix class
//--------------------------------------------------------------------

/**\brief Class representing fixed size symmetric NxN matrix

\tparam T Data type of components
\tparam N Dimensionality
**/
template<class T, int N> 
class sm_matrix_sym {
  static_assert(N>0, "matrix dimension must be >0");
  
  using me_t = sm_matrix_sym<T, N>;
  enum {FLAT_SIZE=(N*(N+1)/2)};

  static int chkidx(int i) 
  {
    assert((i>=0) && (i<N));
    return i;
  }
  static constexpr int index(int i,int j) 
  {
    return  (j <= i) ? (j + (i*(i+1))/2) : (i + (j*(j+1))/2);
  }

  sm_vector<T, FLAT_SIZE> c;
  
  template<class, int> friend class sm_matrix_sym;

  public:
  using elem_t = T;

  sm_matrix_sym()                   noexcept = default;
  sm_matrix_sym(me_t&& that)        noexcept = default;
  sm_matrix_sym(const me_t& that)   noexcept = default;
  sm_matrix_sym& operator=(const me_t& that) noexcept = default;
  sm_matrix_sym& operator=(me_t&& that)      noexcept = default;
  
  template<class T2>
  explicit sm_matrix_sym(const sm_matrix_sym<T2,N>& that) 
  : c{that.c} {}
  
  template<class... An, 
           class Na=detail::fix_num_args<FLAT_SIZE,An...> > 
  explicit sm_matrix_sym(An&&... an) 
  : c{std::forward<An>(an)...} {}

  
  sm_matrix_sym(zero_literal) {c.zero();}
  sm_matrix_sym(one_literal) {diag(1.0);}

  T& operator()(int i,int j) 
  {
    return c(index(chkidx(i),chkidx(j)));
  }
  const T& operator()(int i,int j) const 
  {
    return c(index(chkidx(i),chkidx(j)));
  }

  
  void negate() {
    c.negate();
  }

  template<class T2>
  void operator+=(const sm_matrix_sym<T2,N>& a) {c += a.c;}
  
  template<class T2>
  void operator-=(const sm_matrix_sym<T2,N>& a) {c -= a.c;}
  
  template<class T2>
  void operator*=(T2 z) {c *= z;}
  
  template<class T2>
  void operator/=(T2 z) {c /= z;}

  me_t operator+(const me_t& a) const {
    me_t e{*this};
    e += a; 
    return e;
  }
  me_t operator-(const me_t& a) const {
    me_t e{*this};
    e -= a; 
    return e;
  }
  
  template<class T2>
  me_t operator*(T z) const {
    me_t e{*this}; 
    e *= z;
    return e;
  }
  
  template<class T2>
  me_t operator/(T z) const {
    me_t e{*this};
    e /= z; 
    return e;
  }
  
  me_t operator-() const {
    me_t e{*this};
    e.negate();
    return e;
  }
  
  void diag(T d);

  T bilinear(const sm_vector<T, N>& v, const sm_vector<T, N>& w) const;
  T bilinear(const sm_vector<T, N>& v) const;
};

template<class T, int N>
void sm_matrix_sym<T, N>::diag(T d) {
  for (int i=0;i<N;i++) {
    (*this)(i,i)=d;
    for (int j=0;j<i;j++) {
      (*this)(i,j)=0.0;
    }
  }
}


template<class T, int N>
auto operator*(T a, const sm_matrix_sym<T, N>& m)
-> sm_matrix_sym<T, N>
{
  return m * a;
}

//--------------------------------------------------------------------
//  Raw square matrix class
//--------------------------------------------------------------------

/**\brief Class representing fixed size NxN square matrix

\tparam T Data type of components
\tparam N Dimensionality
**/
template<class T, int N> 
class sm_matrix_sqr {
  static_assert(N>0, "matrix dimension must be >0");
  
  using me_t = sm_matrix_sqr<T, N>;
  enum {FLAT_SIZE = N*N};

  static int chkidx(int i) 
  {
    assert((i>=0) && (i<N));
    return i;
  }
  static constexpr int index(int i,int j) 
  {
    return  N*i+j;
  }
  
  sm_vector<T, FLAT_SIZE> c;
  
  template<class, int> friend class sm_matrix_sqr;


  public:
  using elem_t = T;

  sm_matrix_sqr()                   noexcept = default;
  sm_matrix_sqr(me_t&& that)        noexcept = default;
  sm_matrix_sqr(const me_t& that)   noexcept = default;
  sm_matrix_sqr& operator=(const me_t& that) noexcept = default;
  sm_matrix_sqr& operator=(me_t&& that)      noexcept = default;
  
  
  template<class T2>
  explicit   sm_matrix_sqr(const sm_matrix_sqr<T2,N>& that) 
  : c{that.c} {}
  
  
  template<class... An, 
           class Na=detail::fix_num_args<FLAT_SIZE,An...> > 
  explicit sm_matrix_sqr(An&&... an) 
  : c{std::forward<An>(an)...} {}


  sm_matrix_sqr(zero_literal) {c.zero();}
  sm_matrix_sqr(one_literal) {diag(1.0);}

  T& operator()(int i,int j) 
  {
    return c(index(chkidx(i),chkidx(j)));
  }
  const T& operator()(int i,int j) const 
  {
    return c(index(chkidx(i),chkidx(j)));
  }


  ///Replace matrix by its negative
  void negate() {
    c.negate();
  }

  template<class T2>
  void operator+=(const sm_matrix_sqr<T2,N>& a) {c += a.c;}
  
  template<class T2>
  void operator-=(const sm_matrix_sqr<T2,N>& a) {c -= a.c;}
  
  template<class T2>
  void operator*=(T2 z) {c *= z;}
  
  template<class T2>
  void operator/=(T2 z) {c /= z;}


  me_t operator+(const me_t& a) const {
    me_t e{*this};
    e += a; 
    return e;
  }
  me_t operator-(const me_t& a) const {
    me_t e{*this};
    e -= a; 
    return e;
  }
  
  template<class T2>
  me_t operator*(T2 z) const {
    me_t e{*this}; 
    e *= z;
    return e;
  }
  
  template<class T2>
  me_t operator/(T2 z) const {
    me_t e{*this};
    e /= z; 
    return e;
  }
  
  me_t operator-() const {
    me_t e{*this};
    e.negate();
    return e;
  }

  void diag(T d);

  T bilinear(const sm_vector<T, N>& v, const sm_vector<T, N>& w) const;
  T bilinear(const sm_vector<T, N>& v) const;
};

template<class T, int N>
void sm_matrix_sqr<T, N>::diag(T d) {
  for (int i=0;i<N;i++) {
    (*this)(i,i)=d;
    for (int j=0;j<i;j++) {
      (*this)(i,j)=0.0;
      (*this)(j,i)=0.0;
    }
  }
}

template<class T, int N>
auto operator*(T a, const sm_matrix_sqr<T, N>& m) 
-> sm_matrix_sqr<T, N>
{
  return m * a;
}

//-------------------------------------------------------------------
// Vector times symmetric matrix
//-------------------------------------------------------------------

template<class T, int N>
void sm_vector<T, N>::assign_prod(const sm_matrix_sym<T, N>& m, 
                                  const me_t& w)
{
  for (int i=0; i<N; i++) {
    v[i] = w(0)*m(i,0);
    for (int j=1; j<N; j++) {
      v[i] += w(j) * m(i,j);
    }
  }
}

template<class T, int N>
sm_vector<T, N> operator*(const sm_matrix_sym<T, N>& m, 
                          const sm_vector<T, N>& w)
{
  sm_vector<T, N> erg;
  erg.assign_prod(m, w);
  return erg;
}

template<class T, int N>
sm_vector<T, N> operator*(const sm_vector<T, N>& w, 
                          const sm_matrix_sym<T, N>& m)
{
  sm_vector<T, N> erg;
  erg.assign_prod(w, m);
  return erg;
}

//-----------------------------------------------------------------
// Vector times square matrix
//-----------------------------------------------------------------

template<class T, int N>
void sm_vector<T, N>::assign_prod(const sm_matrix_sqr<T, N>& m, 
                                  const me_t& w)
{
  for (int i=0; i<N; i++) {
    v[i] = w(0)*m(i,0);
    for (int j=1; j<N; j++) {
      v[i] += w(j) * m(i,j);
    }
  }
}

template<class T, int N>
void sm_vector<T, N>::assign_prod(const me_t& w, 
                                  const sm_matrix_sqr<T, N>& m)
{
  for (int i=0; i<N; i++) {
    v[i] = w(0)*m(0,i);
    for (int j=1; j<N; j++) {
      v[i] += w(j) * m(j,i);
    }
  }
}

template<class T, int N>
sm_vector<T, N> operator*(const sm_matrix_sqr<T, N>& m, 
                          const sm_vector<T, N>& w)
{
  sm_vector<T, N> erg;
  erg.assign_prod(m, w);
  return erg;
}

template<class T, int N>
sm_vector<T, N> operator*(const sm_vector<T, N>& w,  
                          const sm_matrix_sqr<T, N>& m)
{
  sm_vector<T, N> erg;
  erg.assign_prod(w, m);
  return erg;
}

//-----------------------------------------------------------------
// Sym matrix as bilinear form 
//-----------------------------------------------------------------

template<class T, int N>
T sm_matrix_sym<T, N>::bilinear(const sm_vector<T, N>& v, 
                                const sm_vector<T, N>& w) const
{
  sm_vector<T, N> t{};   //TODO: exploit symmetry
  t.assign_prod(*this, w);
  return v.dot(t);
}

template<class T, int N>
T sm_matrix_sym<T, N>::bilinear(const sm_vector<T, N>& v) const
{
  const me_t& m=*this;
  T erg = v(0) * v(0) * m(0,0);
  for (int i=1; i<N; i++) {
    T t = m(i,0)*v(0);
    for (int j=1; j<i; j++) {
      t += m(i,j)*v(j);
    }
    erg += v(i) * (v(i)*m(i,i) + 2*t);
  }
  return erg;
}

//-----------------------------------------------------------------
// Square matrix as bilinear form 
//-----------------------------------------------------------------

template<class T, int N>
T sm_matrix_sqr<T, N>::bilinear(const sm_vector<T, N>& v, 
                                const sm_vector<T, N>& w) const
{
  sm_vector<T, N> t;   
  t.assign_prod(*this, w);
  return v.dot(t);
}

template<class T, int N>
T sm_matrix_sqr<T, N>::bilinear(const sm_vector<T, N>& v) const
{
  return bilinear(v,v);
}

//-----------------------------------------------------------------
// Determinant of symmetric 3-matrix
//-----------------------------------------------------------------
template<class T>
T determinant(const sm_matrix_sym<T, 3> &m)
{
  T d = m(0,0) * m(1,1) * m(2,2) 
        + 2 * m(0,1) * m(0,2) * m(1,2)
        - m(0,0) * m(1,2) * m(1,2)
        - m(1,1) * m(0,2) * m(0,2) 
        - m(2,2) * m(0,1) * m(0,1);
  return d;
}

template<class T>
void invert_matrix(const sm_matrix_sym<T, 3>&m, 
                   sm_matrix_sym<T, 3>& erg, T& det)
{
  det = determinant(m);

  erg(0,0) = (-m(1,2)*m(1,2) + m(1,1)*m(2,2) );
  erg(0,1) = ( m(0,2)*m(1,2) - m(0,1)*m(2,2) );
  erg(1,1) = (-m(0,2)*m(0,2) + m(0,0)*m(2,2) );
  erg(0,2) = (-m(0,2)*m(1,1) + m(0,1)*m(1,2) );
  erg(1,2) = ( m(0,1)*m(0,2) - m(0,0)*m(1,2) );
  erg(2,2) = (-m(0,1)*m(0,1) + m(0,0)*m(1,1) );

  erg /= det;
}

//------------------------------------------------------------------
// Determinant of square 2-matrix
//------------------------------------------------------------------
template<class T>
T determinant(const sm_matrix_sqr<T, 2> &m)
{
  T d = m(0,0) * m(1,1) - m(1,0) * m(0,1);
  return d;
}

template<class T>
void invert_matrix(const sm_matrix_sqr<T, 2>&m, 
                   sm_matrix_sqr<T, 2>& erg, T& det)
{
  det = determinant(m);
  erg(0,0) = m(1,1) / det;
  erg(1,1) = m(0,0) / det ;
  erg(0,1) = -m(0,1) / det;
  erg(1,0) = -m(1,0) / det;
}

//------------------------------------------------------------------
//  co/contra-variant vector
//------------------------------------------------------------------

/**\brief Class representing rank-1 tensor

\tparam T Underlying scalar data type
\tparam N Dimensionality
\tparam UP If tensor has upper (true) or lower (false) indices 
**/
template<class T, int N, bool UP> class sm_tensor1 {
  using me_t = sm_tensor1<T, N, UP>;
  using vec_t = sm_vector<T, N>;
  template<class T2> using sim_t = sm_tensor1<T2, N, UP>;

  vec_t c;
  
  template<class, int, bool> friend class sm_tensor1;
  
  public:
  
  using elem_t = T;
  enum {SIZE=N};

  ///Default constructor leaves elements uninitialized.
  sm_tensor1()                      noexcept = default;
  
  sm_tensor1(me_t&& that)           noexcept = default;
  
  ///Copy from rank-1 tensor of same dimension and data type
  sm_tensor1(const me_t& that)      noexcept = default;
  
  ///Assign from rank-1 tensor of same dimension and data type
  sm_tensor1& operator=(const me_t& that) noexcept = default;
  
  sm_tensor1& operator=(me_t&& that)      noexcept = default;
  
  
  ///Copy from tensor with compatible data type
  template<class T2>
  explicit   
  sm_tensor1(const sim_t<T2>& that) noexcept
  : c{that.c} {}
  
  ///Create from single components.
  template<class... An, 
           class Na=detail::fix_num_args<N,An...> > 
  sm_tensor1(An&&... an) 
  : c(std::forward<An>(an)...) {}

  ///Create zero tensor
  sm_tensor1(zero_literal) {
    c.zero();
  }

  ///Access i-th component (counting starts at 0).
  T& operator()(int j) {return c(j);}
  
  ///Get i-th component (counting starts at 0)
  const T& operator()(int j) const {return c(j);}
  
  ///Access underlying raw vector
  vec_t& as_vector() {return c;}
  
  ///Get underlying raw vector
  const vec_t& as_vector() const {return c;}

    
  template<bool UP1>
  void assign_prod(const sm_tensor1<T, N, !UP1>& v, 
                   const sm_tensor2_sym<T, N, UP1, UP>& m) 
  {
    c.assign_prod(v.as_vector(), m.as_matrix());
  }
  
  template<bool UP1>
  void assign_prod(const sm_tensor2_sym<T, N, UP, UP1>& m, 
                   const sm_tensor1<T, N, !UP1>& v) 
  {
    c.assign_prod(m.as_matrix(), v.as_vector());
  }

  ///Replace tensor by its negative
  void negate() {
    c.negate();
  }

  ///Add other compatible tensor
  template<class T2>
  void operator+=(const sim_t<T2>& a) {c += a.c;}
  
  ///Subtract other compatible tensor
  template<class T2>
  void operator-=(const sim_t<T2>& a) {c -= a.c;}
  
  ///Multiply with scalar 
  template<class T2>
  void operator*=(T2 z) {c *= z;}
  
  ///Divide by scalar 
  template<class T2>
  void operator/=(T2 z) {c /= z;}

  ///Compute sum of two tensors
  me_t operator+(const me_t &a) const {
    me_t e{*this};
    e += a; 
    return e;
  }
  
  ///Compute difference of two tensors
  me_t operator-(const me_t &a) const {
    me_t e{*this};
    e -= a; 
    return e;
  }

  ///Compute negative of tensor
  me_t operator-() const {
    me_t e{*this};
    e.negate(); 
    return e;
  }
  
  ///Compute product of tensor with scalar
  me_t operator*(T z) const {
    me_t e{*this};
    e *= z; 
    return e;
  }
  
  ///Compute tensor divided by scalar
  me_t operator/(T z) const {
    me_t e{*this};
    e /= z; 
    return e;
  }
};


//------------------------------------------------------------------
// Scalar * vector
//------------------------------------------------------------------

template<class T, int N, bool UP>
auto inline operator*(T z, const sm_tensor1<T, N, UP>& a) 
->sm_tensor1<T, N, UP>
{
  return a * z;
}

//------------------------------------------------------------------
// Vector-Vector contraction v^i w_i  resp. v_i w^i
//------------------------------------------------------------------
template<class T, int N, bool UP>
T contract(const sm_tensor1<T, N, UP>& v, 
           const sm_tensor1<T, N, !UP>& w)
{
  return v.as_vector().dot(w.as_vector());
} 

template<class T, int N, bool UP>
T operator*(const sm_tensor1<T, N, UP> &v, 
            const sm_tensor1<T, N, !UP> &w)
{
  return contract(v, w);
}

//-----------------------------------------------------------------
// Cross-product for 3-vectors
//-----------------------------------------------------------------

template<class T>
auto cross_product(const sm_tensor1<T,3,true>& a, 
        const sm_tensor1<T,3,true>& b, T vol_elem) 
-> sm_tensor1<T, 3, false>
{
  return {vol_elem * (a(1)*b(2) - a(2)*b(1)),
          vol_elem * (a(2)*b(0) - a(0)*b(2)),
          vol_elem * (a(0)*b(1) - a(1)*b(0))};
}


//-----------------------------------------------------------------
// Symmetric tensors m^ij , m_ij, m^i_j, m_i^j
//-----------------------------------------------------------------

/**\brief Class representing symmetric rank-2 tensor

\tparam T Data type of components
\tparam N Dimensionality
\tparam UP1 If first index is upper (true) or lower (false)
\tparam UP2 If second index is upper (true) or lower (false)
**/
template<class T, int N, bool UP1, bool UP2> 
class sm_tensor2_sym {
  using me_t = sm_tensor2_sym<T, N, UP1, UP2>;
  using mat_t = sm_matrix_sym<T,N>;
  template<class T2> using sim_t = sm_tensor2_sym<T2, N, UP1, UP2>;

  enum {FLAT_SIZE=(N*(N+1)/2)};

  mat_t m;
  
  template<class, int, bool, bool> friend class sm_tensor_sym;
  
  public:
  
  ///Default contructor leaves components uninitialized
  sm_tensor2_sym()                  noexcept = default;
  
  sm_tensor2_sym(me_t&& that)       noexcept = default;
  
  ///Copy from rank-2 tensor of same dimension and data type
  sm_tensor2_sym(const me_t& that)  noexcept = default;
  
  ///Assignment from rank-2 tensor of same dimension and data type
  sm_tensor2_sym& operator=(const me_t& that) noexcept = default;
  
  sm_tensor2_sym& operator=(me_t&& that)      noexcept = default;

  ///Copy from rank-2 tensor of same dimension and compatible data type
  template<class T2> 
  explicit   
  sm_tensor2_sym(const sim_t<T2>& that) 
  : m{that.m} {}

  /** \brief Construct from single components
  
  The order is \f$00, 10,11, \ldots , m0, \dots , mm \f$, where 
  \f$ m=N-1 \f$.
  **/
  template<class... An, 
           class Na=detail::fix_num_args<FLAT_SIZE,An...> > 
  explicit sm_tensor2_sym(An&&... an) 
  : m{std::forward<An>(an)...} {}
  
  
  ///Set to zero
  sm_tensor2_sym(zero_literal) {m.zero();}
  
  ///Set to Kronecker delta
  sm_tensor2_sym(one_literal) {m.diag(1.0);}

  ///Access component \f$ (i,j) \f$ (counting starts at 0)
  T& operator()(int i, int j) {return m(i,j);}
  
  ///Get component \f$ (i,j) \f$ (counting starts at 0)
  const T& operator()(int i, int j) const {return m(i,j);}
  
  ///Access underlying raw matrix of components
  mat_t& as_matrix() {return m;}
  
  ///Get underlying raw matrix of components
  const mat_t& as_matrix() const {return m;}

  ///Replace tensor by its negative
  void negate() {
    m.negate();
  }

  ///Add other compatible tensor
  template<class T2>
  void operator+=(const sim_t<T2>& a) {m += a.m;}
  
  ///Subtract other compatible tensor
  template<class T2>
  void operator-=(const sim_t<T2>& a) {m -= a.m;}
  
  ///Multiply with scalar 
  template<class T2>
  void operator*=(T2 z) {m *= z;}
  
  ///Divide by scalar 
  template<class T2>
  void operator/=(T2 z) {m /= z;}

  ///Compute sum of two tensors
  me_t operator+(const me_t &a) const {
    me_t e{*this};
    e += a; 
    return e;
  }
  
  ///Compute difference of two tensors
  me_t operator-(const me_t &a) const {
    me_t e{*this};
    e -= a; 
    return e;
  }

  ///Compute negative of tensor
  me_t operator-() const {
    me_t e{*this};
    e.negate(); 
    return e;
  }
  
  ///Compute product of tensor with scalar
  me_t operator*(T z) const {
    me_t e{*this};
    e *= z; 
    return e;
  }
  
  ///Compute tensor divided by scalar
  me_t operator/(T z) const {
    me_t e{*this};
    e /= z; 
    return e;
  }

  ///Set component matrix to diagonal matrix
  void diag(const T& d) {
    m.diag(d);
  }
  
  /**\brief Compute contraction with two rank-1 tensors.
  
  @param v Rank-1 tensor to contract first index with
  @param w Rank-1 tensor to contract second index with
  **/
  T contract(const sm_tensor1<T, N, !UP1>& v, 
             const sm_tensor1<T, N, !UP2>& w) const 
  {
    return m.bilinear(v.as_vector(), w.as_vector());
  }
  
  
  /// Contract twice with same rank-1 tensor.
  T quadratic(const sm_tensor1<T, N, !UP1>& v) const 
  {
    return contract_quadratic(*this, v);
  }
};



//-----------------------------------------------------------------
// Scalar * tensor
//-----------------------------------------------------------------

template<class T, int N, bool UP1, bool UP2>
auto operator*(T a, const sm_tensor2_sym<T, N, UP1, UP2> &m) 
-> sm_tensor2_sym<T, N, UP1, UP2> 
{
  return m * a;
}


//-----------------------------------------------------------------
// Contraction m^ij w_j or m_ij w^j or m^i_j w^j
//----------------------------------------------------------------- 


template<class T, int N, bool UP1, bool UP2>
sm_tensor1<T, N, UP1> 
operator*(const sm_tensor2_sym<T, N, UP1, UP2>& m, 
          const sm_tensor1<T, N, !UP2>& w)
{
  sm_tensor1<T, N, UP1> erg;
  erg.assign_prod(m, w);
  return erg;
}

template<class T, int N, bool UP1, bool UP2>
sm_tensor1<T, N, UP2> 
operator*(const sm_tensor1<T, N, !UP1>& w, 
          const sm_tensor2_sym<T, N, UP1, UP2>& m)
{
  sm_tensor1<T, N, UP2> erg;
  erg.assign_prod(w, m);
  return erg;
}


//------------------------------------------------------------------
// Contraction v^i m_ij v^j
//------------------------------------------------------------------


template<class T, int N, bool UP>
T contract_quadratic(const sm_tensor2_sym<T, N, UP, UP>&m, 
                     const sm_tensor1<T, N, !UP> &v) 
{
  return m.as_matrix().bilinear(v.as_vector());
}

//------------------------------------------------------------------
// Determinant
//------------------------------------------------------------------

template<class T, int N, bool UP1, bool UP2>
T determinant(const sm_tensor2_sym<T, N, UP1, UP2> &m)
{
  return determinant(m.as_matrix());
}

//------------------------------------------------------------------
// Metric
//------------------------------------------------------------------

/**\brief Class representing a metric

This allows raising, lowering, and contraction of vectors.
Consists of upper- and lower-index metric tensors and the metric 
determinant.

\tparam T Data type of metric components
\tparam N Dimensionality
**/
template<class T, int N> class sm_metric  {
  public:
  using me_t = sm_metric<T,N>;
  using lo_t = sm_tensor2_sym<T, N, false, false>;
  using up_t = sm_tensor2_sym<T, N, true, true>;
  
  lo_t lo;  
  up_t up;
  T vol_elem, det;

  sm_metric()                       noexcept = default;
  sm_metric(me_t&& that)            noexcept = default;
  sm_metric(const me_t& that)       noexcept = default;
  sm_metric& operator=(const me_t& that) noexcept = default;
  sm_metric& operator=(me_t&& that)      noexcept = default;

  
  ///Copy from metric with compatible component data type
  template<class T2>
  explicit   sm_metric(const sm_metric<T2,N>& that) 
  : lo{that.lo}, up{that.up}, vol_elem{that.vol_elem}, det{that.det} {}

  /**\brief Construct metric
  @param lo_ Lower-index metric tensor      
  @param up_ Upper-index metric tensor      
  @param det_ Determinent of lower-index metric tensor      
  **/
    sm_metric(lo_t lo_, up_t up_, T det_)
  : lo{lo_}, up{up_}, vol_elem{std::sqrt(det_)}, det{det_} {} 

  /**\brief Construct metric
  
  Upper-index metric tensor and determinant will be computed from
  lower-index metric tensor. 
  
  \note Only implemented for N=2,3

  @param lo_ Lower-index metric tensor      
  **/
  explicit sm_metric(lo_t lo_) : lo{lo_} 
  {
    invert_matrix(lo.as_matrix(), up.as_matrix(), det);
    vol_elem = std::sqrt(det);  
  } 

  /**\brief Raise index of vector
  
  @param erg Reference to upper-index vector for storing result
  @param v   Lower-index vector to be raised.
  **/
  void raise(sm_tensor1<T, N, true>& erg, 
             const sm_tensor1<T, N, false>& v) const 
  {
    erg.assign_prod(up, v);
  }

  /**\brief Raise index of vector
  
  @return   Upper-index vector
  @param v  Lower-index vector to be raised.
  **/  
  sm_tensor1<T, N, true> 
  raise(const sm_tensor1<T, N, false>& v) const 
  {
    sm_tensor1<T, N, true> erg{};
    erg.assign_prod(up, v);
    return erg;
  }
  
  
  /**\brief Lower index of vector
  
  @param erg Reference to lower-index vector for storing result
  @param v   Upper-index vector to be lowered.
  **/
  void lower(sm_tensor1<T, N, false>& erg, 
             const sm_tensor1<T, N, true>& v) const 
  {
    erg.assign_prod(lo, v);
  }
  
  /**\brief Lower index of vector
  
  @return   Lower-index vector
  @param v  Upper-index vector to be lowered.
  **/  
  sm_tensor1<T, N, false> lower(const sm_tensor1<T, N, true>& v) const 
  {
    sm_tensor1<T, N, false> erg{};
    erg.assign_prod(lo, v);
    return erg;
  }

  ///Contract two upper-index vectors with metric
  T contract(const sm_tensor1<T, N, true>& v, 
             const sm_tensor1<T, N, true>& w) const 
  {
    return lo.contract(v,w);
  }
  
  ///Contract two lower-index vectors with metric
  T contract(const sm_tensor1<T, N, false>& v, 
             const sm_tensor1<T, N, false>& w) const 
  {
    return up.contract(v,w);
  }
  
  ///Contract upper-index vector with itself using metric 
  T norm2(const sm_tensor1<T, N, true>& v) const {
    return lo.quadratic(v);
  } 
  
  ///Contract lower-index vector with itself using metric
  T norm2(const sm_tensor1<T, N, false>& v) const {
    return up.quadratic(v);
  } 
  
  ///Compute vector norm defined by metric  
  template<bool UP> T norm(const sm_tensor1<T, N, UP>& v) const 
  {
    return std::sqrt(norm2(v));
  }

  /**\brief Compute cross-product of 3-vectors
  
  This is ony defined if N=3, i.e. 3-metric.
  
  @return \f$ \epsilon_{ijk} a^j b^k \f$ where \f$ \epsilon \f$
          is the Levi-Civita tensor
  **/
  auto cross_product(const sm_tensor1<T, 3, true>& a, 
                     const sm_tensor1<T, 3, true>& b) const 
  -> sm_tensor1<T, 3, false>
  {
    static_assert(N==3,"Cross product only defined in 3 dimensions");
    return ::EOS_Toolkit::cross_product(a,b,vol_elem);
  }

  ///Set to spatial part of Minkowski metric
  void minkowski() 
  {
    lo = ONE;
    up = ONE;
    vol_elem = 1.0;
    det      = 1.0;
  }
};


///3-vector with upper indices
using sm_vec3u   = sm_tensor1<real_t, 3, true>;

///3-vector with lower indices
using sm_vec3l   = sm_tensor1<real_t, 3, false>;

///3-dimensional symmetric rank-2 tensor with upper indices
using sm_symt3u  = sm_tensor2_sym<real_t, 3, true, true>;

///3-dimensional symmetric rank-2 tensor with lower indices
using sm_symt3l  = sm_tensor2_sym<real_t, 3, false, false>;

///3-dimensional positive definite metric
using sm_metric3 = sm_metric<real_t, 3>;


}

#endif

