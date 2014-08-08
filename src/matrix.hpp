#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

#include "alignedalloc.hpp"
#include "tri.hpp"

template<typename value_t> 
struct MatrixBlock {
  uint32_t r, c, rows;
  value_t *value;
  uint32_t *col_index;
  uint32_t *cols;

  MatrixBlock() : value(0) {}

  void alloc(uint32_t rows_, uint32_t r_, uint32_t c_,
             size_t col_indexsize, size_t colssize, size_t valuessize) {
    rows = rows_;
    r = r_;
    c = c_;
    AlignmentAllocator<char, 32> alloc;
    value = (value_t *) alloc.allocate(valuessize*sizeof(value_t) +
                                       colssize*sizeof(uint32_t) + 
                                       col_indexsize*sizeof(uint32_t));
    col_index = (uint32_t *) (value + valuessize);
    cols = col_index + col_indexsize;
  }

  ~MatrixBlock() {
    delete [] value;
  }
};

template<typename value_t> 
struct BlockedMatrix {
  MatrixBlock<value_t> *blocks;
  size_t num_blocks;

  BlockedMatrix() : blocks(0) {}

  void alloc(uint32_t num_blocks_) {
    num_blocks = num_blocks_;
    blocks = new MatrixBlock<value_t>[num_blocks];
  }

  ~BlockedMatrix() {
    delete [] blocks;
  }
  MatrixBlock<value_t> &operator[](size_t i) { return blocks[i]; }
  const MatrixBlock<value_t> &operator[](size_t i) const { return blocks[i]; }

  size_t size() const { return num_blocks; }
};

#endif // MATRIX_HPP_INCLUDED
