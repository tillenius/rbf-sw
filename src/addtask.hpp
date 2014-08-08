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
  std::string task_name;

  AdditionTask(const char *name,
               VectorBlockHandle<vec4> &out_,
               const VectorBlockHandle<vec4> &H_,
               const double s_, 
               const VectorBlockHandle<vec4> &D_,
               const bool prio)
   : out(out_), H(H_), s(s_), D(D_)
  {
    register_access(ReadWriteAdd::write, out.handle);
    register_access(ReadWriteAdd::read, H.handle);
    register_access(ReadWriteAdd::read, D.handle);
    task_name = name;
#ifdef USE_MPI
    is_prioritized = prio;
    task_name += (is_prioritized ? "P add" : "add");
#else
    task_name += "add";
#endif
  }
  void run(TaskExecutor<Options> &) {
    for (uint32_t i = 0; i < out.size(); ++i)
      out[i].store_sum_with_scaled_vec(H[i], s, D[i]);
  }
  std::string get_name() { return task_name; }
};

#endif // ADDTASK_HPP_INCLUDED
