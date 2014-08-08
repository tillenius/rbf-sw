#ifndef SIZE_VARS_HPP_INCLUDED
#define SIZE_VARS_HPP_INCLUDED

#include <vector>

// part of variables, separated so that it can be used from
// blocked_vec without cyclic dependencies.

struct size_vars {
  uint32_t nd;
  uint32_t chunk_size;
  uint32_t num_chunks;

  std::vector<size_t> rowrank;
  std::vector<size_t> rowcore;

  uint32_t begin(size_t i) {
    return i*chunk_size;
  }

  uint32_t end(size_t i)  {
    if (i == num_chunks-1)
      return nd;
    return (i+1)*chunk_size;
  }

  size_t chunksize(size_t i) {
    if (i == num_chunks - 1)
      return nd - i*chunk_size;
    return chunk_size;
  }
};

#endif // SIZE_VARS_HPP_INCLUDED
