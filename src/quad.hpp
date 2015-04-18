#ifndef QUAD_HPP_INCLUDED
#define QUAD_HPP_INCLUDED

#include <vector>
#include "alignedalloc.hpp"

#ifdef __linux
#include <x86intrin.h>
#endif

#if defined(_MSC_VER)
#define RBFSW_INLINE __forceinline
#define RBFSW_ALIGN32 __declspec(align(32))
#define RBFSW_ALIGN16 __declspec(align(16))
#else
#define RBFSW_INLINE __attribute__((always_inline))
#define RBFSW_ALIGN32 __attribute__ ((aligned (32)))
#define RBFSW_ALIGN16 __attribute__ ((aligned (16)))
#endif

template <typename value_t>
struct quad {
  value_t data[4];
};


#if defined(__AVX__)
struct vec4 {
  __m256d mdata;

  RBFSW_INLINE const double *elem(size_t i) const { return ((double *) &mdata)+i; }
  RBFSW_INLINE double *elem(size_t i) { return ((double *) &mdata)+i; }

  vec4() {}

  RBFSW_INLINE vec4(double d) {
    mdata = _mm256_set1_pd(d);
  }

  RBFSW_INLINE void load(double *buff) {
    mdata = _mm256_load_pd(buff);
  }

  RBFSW_INLINE void add_scaled_vec(double s, const vec4 &rhs) {
    mdata = _mm256_add_pd(mdata,
                          _mm256_mul_pd(_mm256_set1_pd(s),
                                        rhs.mdata ));
  }
  RBFSW_INLINE void sum_and_store(const vec4 &tmp) {
    mdata = _mm256_add_pd(mdata, tmp.mdata);
  }
  RBFSW_INLINE void store_sum_with_scaled_vec(const vec4 &h, double s, const vec4 &d) {
    mdata = _mm256_add_pd(h.mdata, 
                         _mm256_mul_pd(_mm256_set1_pd(s), d.mdata));
  }
  RBFSW_INLINE void rk4_step(const vec4 &rhs, double s, const vec4 &d1, const vec4 &d2, const vec4 &d3, const vec4 &d4)
  {

    mdata = _mm256_add_pd( rhs.mdata,
                           _mm256_mul_pd(_mm256_set1_pd(s),
                                         _mm256_add_pd( d1.mdata,
                                                        _mm256_add_pd(_mm256_mul_pd( _mm256_set1_pd(2.0),
                                                                                     _mm256_add_pd(d2.mdata, d3.mdata)),
                                                                      d4.mdata))));
  }
} RBFSW_ALIGN32;

#elif defined(__SSE3__)
struct vec4 {
  __m128d mdata[2];

  RBFSW_INLINE const double *elem(size_t i) const { return ((double *) &mdata)+i; }
  RBFSW_INLINE double *elem(size_t i) { return ((double *) &mdata)+i; }

  vec4() {}

  RBFSW_INLINE vec4(double d) {
    mdata[0] = _mm_set1_pd(d);
    mdata[1] = _mm_set1_pd(d);
  }
  RBFSW_INLINE void load(double *buff) {
    mdata[0] = _mm_load_pd(&buff[0]);
    mdata[1] = _mm_load_pd(&buff[2]);
  }
  RBFSW_INLINE void add_scaled_vec(double s, const vec4 &rhs) {
    mdata[0] = _mm_add_pd(mdata[0], _mm_mul_pd(_mm_set1_pd(s), rhs.mdata[0] ));
    mdata[1] = _mm_add_pd(mdata[1], _mm_mul_pd(_mm_set1_pd(s), rhs.mdata[1] ));
  }
  RBFSW_INLINE void sum_and_store(const vec4 &tmp) {
    mdata[0] = _mm_add_pd(tmp.mdata[0], mdata[0] );
    mdata[1] = _mm_add_pd(tmp.mdata[1], mdata[1] );
  }
  RBFSW_INLINE void store_sum_with_scaled_vec(const vec4 &h, double s, const vec4 &d) {
    mdata[0] = _mm_add_pd( h.mdata[0], _mm_mul_pd(_mm_set1_pd(s), d.mdata[0] ));
    mdata[1] = _mm_add_pd( h.mdata[1], _mm_mul_pd(_mm_set1_pd(s), d.mdata[1] ));
  }
  RBFSW_INLINE void rk4_step(const vec4 &rhs, double s, const vec4 &d1, const vec4 &d2, const vec4 &d3, const vec4 &d4) {
     mdata[0] = _mm_add_pd(rhs.mdata[0],
                           _mm_mul_pd(_mm_set1_pd(s),
                                      _mm_add_pd(d1.mdata[0],
                                                 _mm_add_pd(_mm_mul_pd(_mm_set1_pd(2.0),
                                                                       _mm_add_pd(d2.mdata[0], d3.mdata[0])),
                                                            d4.mdata[0]))));
     mdata[1] = _mm_add_pd(rhs.mdata[1],
                           _mm_mul_pd(_mm_set1_pd(s),
                                      _mm_add_pd(d1.mdata[1],
                                                 _mm_add_pd(_mm_mul_pd(_mm_set1_pd(2.0),
                                                                       _mm_add_pd(d2.mdata[1], d3.mdata[1])),
                                                            d4.mdata[1]))));
  }
} RBFSW_ALIGN16;
 
#else

struct vec4 {
  double mdata[4];

  RBFSW_INLINE const double *elem(size_t i) const { return ((double *)&mdata) + i; }
  RBFSW_INLINE double *elem(size_t i) { return ((double *)&mdata) + i; }

  vec4() {}

  RBFSW_INLINE vec4(double d) {
    mdata[0] = d;
    mdata[1] = d;
    mdata[2] = d;
    mdata[3] = d;
  }
  RBFSW_INLINE void load(double *buff) {
    mdata[0] = buff[0];
    mdata[1] = buff[1];
    mdata[2] = buff[2];
    mdata[3] = buff[3];
  }
  RBFSW_INLINE void add_scaled_vec(double s, const vec4 &rhs) {
    mdata[0] += s * rhs.mdata[0];
    mdata[1] += s * rhs.mdata[1];
    mdata[2] += s * rhs.mdata[2];
    mdata[3] += s * rhs.mdata[3];
  }
  RBFSW_INLINE void sum_and_store(const vec4 &tmp) {
    mdata[0] += tmp.mdata[0];
    mdata[1] += tmp.mdata[1];
    mdata[2] += tmp.mdata[2];
    mdata[3] += tmp.mdata[3];
  }
  RBFSW_INLINE void store_sum_with_scaled_vec(const vec4 &h, double s, const vec4 &d) {
    mdata[0] = h.mdata[0] + s * d.mdata[0];
    mdata[1] = h.mdata[1] + s * d.mdata[1];
    mdata[2] = h.mdata[2] + s * d.mdata[2];
    mdata[3] = h.mdata[3] + s * d.mdata[3];
  }
  RBFSW_INLINE void rk4_step(const vec4 &rhs, double s, const vec4 &d1, const vec4 &d2, const vec4 &d3, const vec4 &d4) {
    mdata[0] = rhs.mdata[0] + s*( d1.mdata[0] + 2.0*(d2.mdata[0] + d3.mdata[0]) + d4.mdata[0]);
    mdata[1] = rhs.mdata[1] + s*( d1.mdata[1] + 2.0*(d2.mdata[1] + d3.mdata[1]) + d4.mdata[1]);
    mdata[2] = rhs.mdata[2] + s*( d1.mdata[2] + 2.0*(d2.mdata[2] + d3.mdata[2]) + d4.mdata[2]);
    mdata[3] = rhs.mdata[3] + s*( d1.mdata[3] + 2.0*(d2.mdata[3] + d3.mdata[3]) + d4.mdata[3]);
  }
};
#endif

typedef std::vector<vec4, AlignmentAllocator<vec4, 32> > quadvector_t;

#endif // QUAD_HPP_INCLUDED
