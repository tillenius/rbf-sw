#ifndef QUAD_HPP_INCLUDED
#define QUAD_HPP_INCLUDED

#include <vector>
#include "alignedalloc.hpp"

#ifdef __linux
#include <x86intrin.h>
#endif

template <typename value_t>
struct quad {
  value_t data[4];
};


#if defined(__AVX__)
struct vec4 {
  __m256d mdata;

  const double *elem(size_t i) const __attribute__((always_inline)) { return ((double *) &mdata)+i; }
  double *elem(size_t i) __attribute__((always_inline))  { return ((double *) &mdata)+i; }

  vec4() {}

  vec4(double d) __attribute__((always_inline)) {
    mdata = _mm256_set1_pd(d);
  }

  void load(double *buff) {
    mdata = _mm256_load_pd(buff);
  }

  void add_scaled_vec(double s, const vec4 &rhs) __attribute__((always_inline)) {
    mdata = _mm256_add_pd(mdata,
                          _mm256_mul_pd(_mm256_set1_pd(s),
                                        rhs.mdata ));
  }
  void sum_and_store(const vec4 &tmp) __attribute__((always_inline)) {
    mdata = _mm256_add_pd(mdata, tmp.mdata);
  }
  void store_sum_with_scaled_vec(const vec4 &h, double s, const vec4 &d) __attribute__((always_inline)) {
    mdata = _mm256_add_pd(h.mdata, 
                         _mm256_mul_pd(_mm256_set1_pd(s), d.mdata));
  }
  void rk4_step(const vec4 &rhs, double s, const vec4 &d1, const vec4 &d2, const vec4 &d3, const vec4 &d4)
  __attribute__((always_inline)) {

    mdata = _mm256_add_pd( rhs.mdata,
                           _mm256_mul_pd(_mm256_set1_pd(s),
                                         _mm256_add_pd( d1.mdata,
                                                        _mm256_add_pd(_mm256_mul_pd( _mm256_set1_pd(2.0),
                                                                                     _mm256_add_pd(d2.mdata, d3.mdata)),
                                                                      d4.mdata))));
  }
} __attribute__ ((aligned (32)));

#elif defined(__SSE3__)
struct vec4 {
  __m128d mdata[2];

  const double *elem(size_t i) const __attribute__((always_inline)) { return ((double *) &mdata)+i; }
  double *elem(size_t i) __attribute__((always_inline))  { return ((double *) &mdata)+i; }

  vec4() {}

  vec4(double d) __attribute__((always_inline)) {
    mdata[0] = _mm_set1_pd(d);
    mdata[1] = _mm_set1_pd(d);
  }
  void load(double *buff) {
    mdata[0] = _mm_load_pd(&buff[0]);
    mdata[1] = _mm_load_pd(&buff[2]);
  }
  void add_scaled_vec(double s, const vec4 &rhs) __attribute__((always_inline)) {
    mdata[0] = _mm_add_pd(mdata[0], _mm_mul_pd(_mm_set1_pd(s), rhs.mdata[0] ));
    mdata[1] = _mm_add_pd(mdata[1], _mm_mul_pd(_mm_set1_pd(s), rhs.mdata[1] ));
  }
  void sum_and_store(const vec4 &tmp) __attribute__((always_inline)) {
    mdata[0] = _mm_add_pd(tmp.mdata[0], mdata[0] );
    mdata[1] = _mm_add_pd(tmp.mdata[1], mdata[1] );
  }
  void store_sum_with_scaled_vec(const vec4 &h, double s, const vec4 &d) __attribute__((always_inline)) {
    mdata[0] = _mm_add_pd( h.mdata[0], _mm_mul_pd(_mm_set1_pd(s), d.mdata[0] ));
    mdata[1] = _mm_add_pd( h.mdata[1], _mm_mul_pd(_mm_set1_pd(s), d.mdata[1] ));
  }
  void rk4_step(const vec4 &rhs, double s, const vec4 &d1, const vec4 &d2, const vec4 &d3, const vec4 &d4) __attribute__((always_inline)) {
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
} __attribute__ ((aligned (16)));
 
#else

struct vec4 {
  double mdata[4];

  const double *elem(size_t i) const __attribute__((always_inline)) { return ((double *) &mdata)+i; }
  double *elem(size_t i) __attribute__((always_inline))  { return ((double *) &mdata)+i; }

  vec4() {}

  vec4(double d) __attribute__((always_inline)) {
    mdata[0] = d;
    mdata[1] = d;
    mdata[2] = d;
    mdata[3] = d;
  }
  void load(double *buff) {
    mdata[0] = buff[0];
    mdata[1] = buff[1];
    mdata[2] = buff[2];
    mdata[3] = buff[3];
  }
  void add_scaled_vec(double s, const vec4 &rhs) __attribute__((always_inline)) {
    mdata[0] += s * rhs.mdata[0];
    mdata[1] += s * rhs.mdata[1];
    mdata[2] += s * rhs.mdata[2];
    mdata[3] += s * rhs.mdata[3];
  }
  void sum_and_store(const vec4 &tmp) {
    mdata[0] += tmp.mdata[0];
    mdata[1] += tmp.mdata[1];
    mdata[2] += tmp.mdata[2];
    mdata[3] += tmp.mdata[3];
  }
  void store_sum_with_scaled_vec(const vec4 &h, double s, const vec4 &d) __attribute__((always_inline)) {
    mdata[0] = h.mdata[0] + s * d.mdata[0];
    mdata[1] = h.mdata[1] + s * d.mdata[1];
    mdata[2] = h.mdata[2] + s * d.mdata[2];
    mdata[3] = h.mdata[3] + s * d.mdata[3];
  }
  void rk4_step(const vec4 &rhs, double s, const vec4 &d1, const vec4 &d2, const vec4 &d3, const vec4 &d4) __attribute__((always_inline)) {
    mdata[0] = rhs.mdata[0] + s*( d1.mdata[0] + 2.0*(d2.mdata[0] + d3.mdata[0]) + d4.mdata[0]);
    mdata[1] = rhs.mdata[1] + s*( d1.mdata[1] + 2.0*(d2.mdata[1] + d3.mdata[1]) + d4.mdata[1]);
    mdata[2] = rhs.mdata[2] + s*( d1.mdata[2] + 2.0*(d2.mdata[2] + d3.mdata[2]) + d4.mdata[2]);
    mdata[3] = rhs.mdata[3] + s*( d1.mdata[3] + 2.0*(d2.mdata[3] + d3.mdata[3]) + d4.mdata[3]);
  }
};
#endif

typedef std::vector<vec4, AlignmentAllocator<vec4, 32> > quadvector_t;

#endif // QUAD_HPP_INCLUDED
