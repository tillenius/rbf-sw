#ifndef VARIABLES_HPP_INCLUDED
#define VARIABLES_HPP_INCLUDED

#include "blocked_vec.hpp"
#include "size_vars.hpp"

struct atmdata {
  double f;
  double x, y, z;
  double p_u[3], p_v[3], p_w[3];
  double ghm;
  double gradghm[3];
};

struct rowcomp {
  bool operator()(const std::pair<uint32_t, uint32_t> &l, size_t r) { return l.first < r; }
};

// todo: separate class this into data and logic
struct variables : public size_vars {

  double gh0;

  BlockedMatrix< quad<double> > D;
  BlockedVector<atmdata> atm;

private:

  void loadFull(const char *name, size_t &m, size_t &n, std::vector<double> &out) {
    // to support large files we should consider memory-mapping the file instead.
    // but that cannot easily be done in a portable way, so I keep loading the whole
    // file for now, until the problem needs to be handled.
    FILE *f=fopen(name, "r");
    if (f == NULL) { std::cerr<<"Error opening file " << name << std::endl; exit(1); }
    fread(&m, sizeof(size_t), 1, f) && 1;
    fread(&n, sizeof(size_t), 1, f) && 1;
    out.resize(m*n);
    fread( &out[0], sizeof(double), m*n, f) && 1;
    fclose(f);
  }

  void loadSparse(const char *name,
                  std::vector< std::pair<uint32_t, uint32_t> > &idx,
                  std::vector< quad<double> > &data)
  {
    // to support large files we should consider memory-mapping the file instead.
    // but that cannot easily be done in a portable way, so I keep loading the whole
    // file for now, until the problem needs to be handled.
    FILE *f = fopen(name, "r");
    if (f == NULL) { std::cerr<<"Error opening file " << name << std::endl; exit(1); }

    size_t mn, nnz;
    fread(&mn, sizeof(size_t), 1, f) && 1;
    fread(&nnz, sizeof(size_t), 1, f) && 1;

    idx.resize(nnz);
    fread(&idx[0], sizeof(std::pair<uint32_t, uint32_t>), nnz, f) && 1;

    data.resize(nnz);
    fread(&data[0], sizeof(quad<double>), nnz, f) && 1;

    fclose(f);
  }

  struct ConvQuadVecTask : public Task<Options> {
    size_t size;
    double *src;
    VectorBlockHandle<vec4> * &dst;
    size_t rank;

    ConvQuadVecTask(size_t size_,
                    double *src_,
                    VectorBlockHandle<vec4> * &dst_,
                    Handle<Options> &h,
                    size_t rank_)
      : size(size_), src(src_), dst(dst_), rank(rank_) {
      this->time = -6.0;
      register_access(ReadWriteAdd::add, h);
    }
    void run(TaskExecutor<Options> &) {
      dst = new VectorBlockHandle<vec4>();
      dst->alloc(size, rank);
      memcpy(dst->get_data(), src, size*sizeof(vec4));
    }
  };

  struct ConvAtmVecTask : public Task<Options> {
    size_t size;
    double *src;
    VectorBlock<atmdata> &dst;

    ConvAtmVecTask(size_t size_,
                   double *src_,
                   VectorBlock<atmdata> &dst_,
                   Handle<Options> &h)
    : size(size_), src(src_), dst(dst_) {
      this->time = -7.0;
      register_access(ReadWriteAdd::add, h);
    }
    void run(TaskExecutor<Options> &) {
      dst.alloc(size);
      memcpy(dst.value, src, size*sizeof(atmdata));
    }
  };

  // ConvInfo: Data and some utility function shared between all S
  struct ConvInfo {
    variables &vars;
    std::vector< std::pair<uint32_t, uint32_t> > idx;
    std::vector< quad<double> > data;

    ConvInfo(variables &vars_) : vars(vars_) {}
  };

  void find_distribution(size_t num_cores, size_t num_ranks) {
    rowrank.resize(num_chunks);
    rowcore.resize(num_chunks);

    for (size_t i = 0; i < num_chunks; ++i) {
      rowrank[i] = (i*num_ranks) / (num_chunks);
      rowcore[i] = (i*num_ranks*num_cores) / (num_chunks) - rowrank[i]*num_cores;
    }
  }

  struct BlockDTask : public Task<Options> {
    size_t brow, bcol, Didx;
    ConvInfo &info;

    BlockDTask(size_t brow_, size_t bcol_, size_t Didx_,
               ConvInfo &info_, 
               Handle<Options> &h)
    : brow(brow_), bcol(bcol_), Didx(Didx_), info(info_)
    {
      this->time = -8.0;
      register_access(ReadWriteAdd::add, h);
      // Note: D[Didx] does not exist at construction time, but only at execution time.
    }

    void run(TaskExecutor<Options> &) {

      const size_t rbegin(info.vars.begin(brow));
      const size_t rend(info.vars.end(brow));
      const size_t cbegin(info.vars.begin(bcol));
      const size_t cend(info.vars.end(bcol));

      std::vector<uint32_t> col_index;
      std::vector<uint32_t> cols;
      std::vector< quad<double> > value;

      uint32_t curr_row = 0;
      for (size_t i = std::lower_bound(info.idx.begin(), info.idx.end(), rbegin, rowcomp()) - info.idx.begin(); i != info.idx.size(); ++i) {
        const uint32_t row = info.idx[i].first;
        if (row >= rend)
          break;

        const uint32_t col = info.idx[i].second;
        if (col < cbegin || col >= cend)
          continue;

        const quad<double> &val(info.data[i]);

        const uint32_t localrow = row - rbegin;
        const uint32_t localcol = col - cbegin;

        for (; curr_row <= localrow; ++curr_row)
          col_index.push_back(cols.size());
        cols.push_back(localcol);
        value.push_back(val);
      }
      for (; curr_row <= rend-rbegin; ++curr_row)
        col_index.push_back(cols.size());

      MatrixBlock< quad<double> > &dst(info.vars.D[Didx]);

      dst.alloc(rend-rbegin, brow, bcol, col_index.size(), cols.size(), value.size()); // cols.size == value.size

      memcpy(dst.col_index, &col_index[0], col_index.size() * sizeof(uint32_t));
      std::vector<uint32_t>().swap(col_index);

      memcpy(dst.cols, &cols[0], cols.size() * sizeof(uint32_t));
      std::vector<uint32_t>().swap(cols);

      memcpy(dst.value, &value[0], value.size() * sizeof(quad<double>));
      std::vector< quad<double> >().swap(value);
    }
  };

  size_t get_num_chunks(size_t nd, size_t chunk_size) {
    size_t num_chunks = nd / chunk_size;
    if (num_chunks == 0)
      num_chunks = 1;
    return num_chunks;
  }

  // load "D"
  void load(tasklib_t &tl, const char *name, BlockedMatrix< quad<double> > &D, size_t rank) {
    ConvInfo info(*this);

    loadSparse(name, info.idx, info.data);

    Handle<Options> h;
    std::vector<Task<Options> *> tasks;

    uint32_t curr_brow = 0;
    std::vector<size_t> elemcount(num_chunks);
    for (size_t i = 0; i < info.idx.size(); ++i) {
      const size_t row = info.idx[i].first;
      const size_t col = info.idx[i].second;
      const size_t brow = std::min<size_t>(row / chunk_size, num_chunks-1);
      const size_t bcol = std::min<size_t>(col / chunk_size, num_chunks-1);

      if (brow != curr_brow) {
        std::vector<size_t>(num_chunks).swap(elemcount);
        curr_brow = brow;
      }

      if (elemcount[bcol] == 0) {
        elemcount[bcol] = 1;
        Task<Options> *task =
            new BlockDTask(brow, bcol, tasks.size(), info, h);
        task->pinned_to = rowcore[brow];
        tasks.push_back(task);
      }
    }

    D.alloc(tasks.size());
    for (size_t i = 0; i < tasks.size(); ++i)
      tl.submit(tasks[i], tasks[i]->pinned_to);
    tl.wait(h);
  }

  // load "atm"
  void load(tasklib_t &tl, const char *name, BlockedVector< atmdata > &out) {
    size_t m, n;
    std::vector<double> data;
    loadFull(name, m, n, data);
    out.alloc(num_chunks);

    Handle<Options> h;

    for (size_t i = 0; i < num_chunks; ++i) {
      const size_t size = chunksize(i);
      Task<Options> *task = new ConvAtmVecTask(size, &data[i*chunk_size*m], out[i], h);
      task->pinned_to = rowcore[i];
      tl.submit(task, rowcore[i]);
    }
    tl.wait(h);
  }

public:

  // load "H"
  void load(tasklib_t &tl, const char *name, BlockedVectorHandle<vec4> &out) {
    size_t m, n;
    std::vector<double> data;
    loadFull(name, m, n, data);
    out.alloc(num_chunks);

    Handle<Options> h;

    for (size_t i = 0; i < num_chunks; ++i) {
      const size_t size = chunksize(i);
      Task<Options> *task = new ConvQuadVecTask(size, &data[i*chunk_size*m], out.chunks[i], h, rowrank[i]);
      task->pinned_to = rowcore[i];
      tl.submit(task, rowcore[i]);
    }
    tl.wait(h);
  }

  void init(tasklib_t &tl, const char *mldatapath_, size_t num_cores, size_t num_ranks, size_t rank, size_t chunk_size_) {

    chunk_size = chunk_size_;

    std::string mldatapath(mldatapath_);

    // load scalars from "params": nd, atm_gh0
    std::string paramname(mldatapath + "/params");
    FILE *f = fopen(paramname.c_str(), "rb");
    if (f == NULL) { std::cerr<<"Error opening '"<<paramname<<"'" << std::endl; exit(1); }
    uint64_t ndtmp;
    assert(fread(&ndtmp, sizeof(uint64_t), 1, f) == 1);
    assert(fread(&gh0, sizeof(double), 1, f) == 1);
    fclose(f);
    nd = (uint32_t) ndtmp;

    // setup class
    num_chunks = get_num_chunks(nd, chunk_size);
    find_distribution(num_cores, num_ranks);

    load(tl, (mldatapath + "/D").c_str(), D, rank);
    load(tl, (mldatapath + "/atm").c_str(), atm);
  }
};

#endif // VARIABLES_HPP_INCLUDED
