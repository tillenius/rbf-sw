#ifndef ADDTASK_HPP_INCLUDED
#define ADDTASK_HPP_INCLUDED

//=============================================================================
// Solver::AdditionTask
//=============================================================================
class AdditionTask : public sgmpi::MPITask<Options> {
public:
  VectorBlockHandle<vec4> &out;
  const VectorBlockHandle<vec4> &H;
  const double s;
  const VectorBlockHandle<vec4> &D;

  AdditionTask(const double time_,
               VectorBlockHandle<vec4> &out_,
               const VectorBlockHandle<vec4> &H_,
               const double s_, 
               const VectorBlockHandle<vec4> &D_,
               const int r)
   : out(out_), H(H_), s(s_), D(D_)
  {
    this->time = time_;
    register_access(ReadWriteAdd::write, out.handle);
    register_access(ReadWriteAdd::read, H.handle);
    register_access(ReadWriteAdd::read, D.handle);
#ifdef USE_MPI
    is_prioritized = prio_row[r] == 1;
#endif
  }
  void run(TaskExecutor<Options> &) {
    for (uint32_t i = 0; i < out.size(); ++i)
      out[i].store_sum_with_scaled_vec(H[i], s, D[i]);
  }
};

#endif // ADDTASK_HPP_INCLUDED
