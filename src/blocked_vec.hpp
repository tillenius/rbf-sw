#ifndef BLOCKED_VEC_HPP_INCLUDED
#define BLOCKED_VEC_HPP_INCLUDED

#include "quad.hpp"
#include "alignedalloc.hpp"
#include "matrix.hpp"
#include "size_vars.hpp"
#include "tasklib.hpp"

// VectorBlock: a chunk in a BlockedVector
template<typename value_t> 
struct VectorBlock {
  uint32_t size;
  value_t *value;

  VectorBlock() : value(0) {}

  void alloc(uint32_t size_) {
    size = size_;
    AlignmentAllocator<value_t, 32> alloc;
    value = alloc.allocate(size);
  }

  ~VectorBlock() {
    delete [] value;
  }
  value_t &operator[](size_t i) { return value[i]; }
  const value_t &operator[](size_t i) const { return value[i]; }
  value_t &operator()(size_t i) { return value[i]; }
  const value_t &operator()(size_t i) const { return value[i]; }
};

// BlockedVector: a vector built from VectorBlocks
template<typename value_t> 
struct BlockedVector {
  VectorBlock<value_t> *blocks;
  size_t num_blocks;

  BlockedVector() : blocks(0) {}

  void alloc(uint32_t num_blocks_) {
    num_blocks = num_blocks_;
    blocks = new VectorBlock<value_t>[num_blocks];
  }

  ~BlockedVector() {
    delete [] blocks;
  }
  VectorBlock<value_t> &operator[](size_t i) { return blocks[i]; }
  const VectorBlock<value_t> &operator[](size_t i) const { return blocks[i]; }
};



// todo: move VectorBlockHandle and BlockedVectorHandle to own file?

// VectorBlockHandle: chunk of a BlockedVectorHandle. Has an associated handle.
template<typename value_t>
struct VectorBlockHandle {
  mutable sgmpi::MPIHandle<Options> handle;

#ifdef USE_MPI
#else
  VectorBlockHandle() : data(0) {}
#endif

  void alloc(uint32_t size_, int rowrank) {
    AlignmentAllocator<value_t, 32> alloc;

#ifdef USE_MPI
    handle.data = (double *) alloc.allocate(size_);
    handle.size = size_*sizeof(value_t)/sizeof(double);
    handle.set_rank(rowrank);
#else
    data = alloc.allocate(size_);
    data_size = size_;
#endif
  }

  double operator()(uint32_t index, uint32_t quadidx) const;
  double x(uint32_t index, uint32_t quadidx) const;
  double y(uint32_t index, uint32_t quadidx) const;
  double z(uint32_t index, uint32_t quadidx) const;
  double l(uint32_t index, uint32_t quadidx) const;

  double &operator()(uint32_t index, uint32_t quadidx);

  value_t &operator[](uint32_t index) { return get_data()[index]; }
  value_t &operator()(uint32_t index) { return get_data()[index]; }
  const value_t &operator[](uint32_t index) const { return get_data()[index]; }
  const value_t &operator()(uint32_t index) const { return get_data()[index]; }

#ifdef USE_MPI
  size_t size() const { return handle.size/(sizeof(value_t)/sizeof(double)); }
  value_t *get_data() const { return (value_t *) handle.data; }
#else
  value_t *data;
  size_t data_size;
  size_t size() const { return data_size; }
  value_t *get_data() const { return data; }
#endif
};

// BlockedVectorHandle: a vector built from VectorBlockHandles. (each chunk has an associated handle)
template<typename value_t>
struct BlockedVectorHandle { 
  typedef typename Options::AccessInfoType::Type accesstype_t;

  std::vector<VectorBlockHandle<value_t> *> chunks;

  void alloc(size_t num_chunks) {
    chunks.resize(num_chunks);
  }

  void init(tasklib_t &tl, size_vars &vars, const char *name) {
    alloc(vars.num_chunks);

    // todo: this part should be moved out somewhere

    Handle<Options> h;
    for (size_t i = 0; i < vars.num_chunks; ++i) {
      const size_t size = vars.chunksize(i);

      struct InitBlockTask : public Task<Options> {
        size_t size;
        VectorBlockHandle<value_t> * &chunk;
        size_t rank;
        InitBlockTask(size_t size_, VectorBlockHandle<value_t> * &chunk_, Handle<Options> &h, size_t rank_)
        : size(size_), chunk(chunk_), rank(rank_) {
          register_access(ReadWriteAdd::read, h);
        }
        void run(TaskExecutor<Options> &) {
          chunk = new VectorBlockHandle<value_t>();
          chunk->alloc(size, rank);
        }
      };

      Task<Options> *task = new InitBlockTask(size, chunks[i], h, vars.rowrank[i]);
      task->pinned_to = vars.rowcore[i];
      tl.submit(task, vars.rowcore[i]);
    }
    tl.wait(h);
  }

  template<typename TaskType>
  void register_access(TaskType &task, accesstype_t at) const {
    for (size_t i = 0; i < chunks.size(); ++i)
      task.register_access(at, chunks[i]->handle);
  }

  template<typename TaskType>
  void register_access(TaskType &task, accesstype_t at, uint32_t index) const {
    if (chunks.size() == 1)
      task.register_access(at, chunks[0]->handle);
    else if (index == ~static_cast<uint32_t>(0)) {
      for (size_t i = 0; i < chunks.size(); ++i)
        task.register_access(at, chunks[i]->handle);
    }
    else
      task.register_access(at, chunks[index]->handle);
  }

#ifdef USE_MPI
  VectorBlockHandle<value_t> &register_access_any(int rank, Task<Options> &task, accesstype_t at, uint32_t index) const {
    for (size_t i = 0; i < chunks.size(); ++i) {
      if (chunks[i]->handle.get_rank() == rank) {
        task.register_access(at, chunks[i]->handle);
        return *chunks[i];
      }
    }
    assert(false);
    return *chunks[0];
  }
#endif

  const VectorBlockHandle<value_t>& operator()(uint32_t index) const { return *chunks[index]; }
  VectorBlockHandle<value_t>& operator()(uint32_t index) { return *chunks[index]; }

  const VectorBlockHandle<value_t>& operator[](uint32_t index) const { return *chunks[index]; }
  VectorBlockHandle<value_t>& operator[](uint32_t index) { return *chunks[index]; }

  size_t size() const { return chunks.size(); }
};


template<> double VectorBlockHandle< vec4 >::operator()(uint32_t index, uint32_t quadidx) const { return *get_data()[index].elem(quadidx); }
template<> double VectorBlockHandle< quad<vec4> >::x(uint32_t index, uint32_t quadidx) const { return *get_data()[index].data[0].elem(quadidx); }
template<> double VectorBlockHandle< quad<vec4> >::y(uint32_t index, uint32_t quadidx) const { return *get_data()[index].data[1].elem(quadidx); }
template<> double VectorBlockHandle< quad<vec4> >::z(uint32_t index, uint32_t quadidx) const { return *get_data()[index].data[2].elem(quadidx); }
template<> double VectorBlockHandle< quad<vec4> >::l(uint32_t index, uint32_t quadidx) const { return *get_data()[index].data[3].elem(quadidx); }
template<> double &VectorBlockHandle< vec4 >::operator()(uint32_t index, uint32_t quadidx) { return *get_data()[index].elem(quadidx); }









template<typename Type>
void dump(BlockedVectorHandle< Type > &T, const char *name) {
  FILE *fout = fopen(name, "wb");
  for (size_t i = 0; i < T.chunks.size(); ++i) {
    for (size_t j = 0; j < T.chunks[i]->size(); ++j) {
      VectorBlockHandle<Type> &b(*T.chunks[i]);
      myfprintf(fout, b.data[j]);
    }
  }
  fclose(fout);
}

#endif // BLOCKED_VEC_HPP_INCLUDED
