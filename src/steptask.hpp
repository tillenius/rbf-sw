#ifndef STEPTASK_HPP_INCLUDED
#define STEPTASK_HPP_INCLUDED

//=============================================================================
// Solver::StepTask
//=============================================================================
class StepTask : public sgmpi::MPITask<Options> {
private:
  VectorBlockHandle<vec4> &Hdst, &Hsrc;
  const double s;
  const VectorBlockHandle<vec4> &d1, &d2, &d3, &d4;

public:

  StepTask(const double time_,
           VectorBlockHandle<vec4> &Hdst_,
           VectorBlockHandle<vec4> &Hsrc_,
           const double s_,
           const VectorBlockHandle<vec4> &d1_,
           const VectorBlockHandle<vec4> &d2_,
           const VectorBlockHandle<vec4> &d3_,
           const VectorBlockHandle<vec4> &d4_,
           const int r)
  : Hdst(Hdst_), Hsrc(Hsrc_), s(s_), d1(d1_), d2(d2_), d3(d3_), d4(d4_)
  {
    this->time = time_;
    register_access(ReadWriteAdd::read, d4.handle);
    register_access(ReadWriteAdd::read, d3.handle);
    register_access(ReadWriteAdd::read, d2.handle);
    register_access(ReadWriteAdd::read, d1.handle);
    if (&Hdst != &Hsrc) {
      register_access(ReadWriteAdd::read, Hsrc.handle);
      register_access(ReadWriteAdd::write, Hdst.handle);
    }
    else {
      register_access(ReadWriteAdd::add, Hdst.handle);
    }

#ifdef USE_MPI
    is_prioritized = prio_row[r] == 1;
#endif
  }

  void run(TaskExecutor<Options> &) {
    for (uint32_t r = 0; r < Hsrc.size(); ++r)
      Hdst[r].rk4_step(Hsrc[r], s, d1[r], d2[r], d3[r], d4[r]);
  }
};

#endif // STEPTASK_HPP_INCLUDED
