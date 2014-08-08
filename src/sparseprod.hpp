#ifndef SPARSEPROD_HPP_INCLUDED
#define SPARSEPROD_HPP_INCLUDED

extern std::vector<uint32_t> prio_row;

//=============================================================================
// sparseProdTask
// T(r) = T(r) + DP(r, c) * H(c)
//=============================================================================
class SparseProdTask : public sgmpi::MPITask<Options> {
public:
  VectorBlockHandle< quad<vec4> > &T;
  const MatrixBlock< quad<double> > &DP;
  const VectorBlockHandle<vec4> &H;
  std::string task_name;
  bool overwrite;

  SparseProdTask(const char *name,
                 uint32_t r, uint32_t c,
                 VectorBlockHandle< quad<vec4> > &T_,
                 const MatrixBlock< quad<double> > &DP_,
                 const VectorBlockHandle<vec4> &H_,
                 bool overwrite_)
  : T(T_), DP(DP_), H(H_), overwrite(overwrite_)
  {
    if (overwrite)
      register_access(ReadWriteAdd::write, T.handle);
    else
      register_access(ReadWriteAdd::add, T.handle);
    register_access(ReadWriteAdd::read, H.handle);
    task_name = name;
#ifdef USE_MPI
    is_prioritized = prio_row[r];
    task_name += (is_prioritized ? "P sparseprod" : "sparseprod");
#else
    task_name += "sparseprod";
#endif
  }

  template<bool overwrite>
  void run_internal() {
    size_t value_index = 0;
    for (uint32_t row = 0; row < DP.rows; ++row) {
      const uint32_t starting_col_index = DP.col_index[row];
      const uint32_t stopping_col_index = DP.col_index[row+1];
      if (starting_col_index == stopping_col_index) {
        if (overwrite) {
          T(row).data[0] = 0.0;
          T(row).data[1] = 0.0;
          T(row).data[2] = 0.0;
          T(row).data[3] = 0.0;
        }
        continue;
      }

      register vec4 tmp[4] = {0.0, 0.0, 0.0, 0.0};

      for (uint32_t cidx = starting_col_index; cidx < stopping_col_index; ++cidx) {
        const uint32_t col = DP.cols[cidx];
        quad<double> dval = DP.value[value_index++];
        tmp[0].add_scaled_vec(dval.data[0], H[col]);
        tmp[1].add_scaled_vec(dval.data[1], H[col]);
        tmp[2].add_scaled_vec(dval.data[2], H[col]);
        tmp[3].add_scaled_vec(dval.data[3], H[col]);
      }

      if (overwrite) {
        T(row).data[0] = tmp[0];
        T(row).data[1] = tmp[1];
        T(row).data[2] = tmp[2];
        T(row).data[3] = tmp[3];
      }
      else {
        T(row).data[0].sum_and_store(tmp[0]);
        T(row).data[1].sum_and_store(tmp[1]);
        T(row).data[2].sum_and_store(tmp[2]);
        T(row).data[3].sum_and_store(tmp[3]);
      }
    }
  }

  void run(TaskExecutor<Options> &) { 
    if (overwrite)
      run_internal<true>();
    else
      run_internal<false>();
  }

  std::string get_name() { return task_name; }
};

#endif // SPARSEPROD_HPP_INCLUDED
