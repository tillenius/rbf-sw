#ifndef RHS_HPP_INCLUDED
#define RHS_HPP_INCLUDED

#include "variables.hpp"
#include "sparseprod.hpp"
#include "tasklib.hpp"

//=============================================================================
// RHS::CalcRHSTask
//=============================================================================
class CalcRHSTask : public sgmpi::MPITask<Options> {
public:
  const VectorBlock<atmdata> &atm;
  const VectorBlockHandle<vec4> &H;
  const VectorBlockHandle< quad<vec4> > &T;
  VectorBlockHandle<vec4> &F;
  const double gh0;
  std::string task_name;

  CalcRHSTask(const char *name,
              const VectorBlock<atmdata> &atm_,
              const VectorBlockHandle<vec4> &H_,
              const VectorBlockHandle< quad<vec4> > &T_,
              const double gh0_,
              VectorBlockHandle<vec4> &F_,
              bool prio)
  : atm(atm_), H(H_), T(T_), F(F_), gh0(gh0_)
  {
    register_access(ReadWriteAdd::write, F.handle);
    register_access(ReadWriteAdd::read, H.handle);
    register_access(ReadWriteAdd::read, T.handle);

    task_name = name;
#ifdef USE_MPI
    is_prioritized = prio;
    task_name += (is_prioritized ? "P rhs" : "rhs" );
#else
    task_name += "rhs";
#endif
  }

  void run(TaskExecutor<Options> &) {

    for (uint32_t i = 0; i < H.size(); ++i) {
      const atmdata &a(atm[i]);

      const double p = -(  H(i,0) * T.x(i,0)
                         + H(i,1) * T.y(i,0)
                         + H(i,2) * T.z(i,0)
                         + a.f * (a.y * H(i,2) - a.z * H(i,1)) + T.x(i,3));
      const double q = -(  H(i,0) * T.x(i,1)
                         + H(i,1) * T.y(i,1)
                         + H(i,2) * T.z(i,1)
                         + a.f * (a.z * H(i,0) - a.x * H(i,2)) + T.y(i,3));
      const double s = -(  H(i,0) * T.x(i,2)
                         + H(i,1) * T.y(i,2)
                         + H(i,2) * T.z(i,2)
                         + a.f * (a.x * H(i,1) - a.y * H(i,0)) + T.z(i,3));

      F(i,0) = a.p_u[0]*p + a.p_u[1]*q + a.p_u[2]*s + T.l(i,0);
      F(i,1) = a.p_v[0]*p + a.p_v[1]*q + a.p_v[2]*s + T.l(i,1);
      F(i,2) = a.p_w[0]*p + a.p_w[1]*q + a.p_w[2]*s + T.l(i,2);

      F(i,3) = -(  H(i,0) * (T.x(i,3) - a.gradghm[0])
                 + H(i,1) * (T.y(i,3) - a.gradghm[1])
                 + H(i,2) * (T.z(i,3) - a.gradghm[2])
                 + (H(i,3)+gh0-a.ghm) * (T.x(i,0) + T.y(i,1) + T.z(i,2)))
               + T.l(i,3);
    }
  }
  std::string get_name() { return task_name; }
};

//=============================================================================
// RHS
//=============================================================================
struct RHS {
  BlockedVectorHandle< quad<vec4> > T;
  variables &vars;

  RHS(variables &vars_) : vars(vars_) {
    T.init(*tl, vars, "T");
  }

  void calc_rhs(const char *name,
                double t,
                BlockedVectorHandle<vec4> &H,
                BlockedVectorHandle<vec4> &F) {

    // T(r) = T(r) + DP(r, c) * H(c) / a
    uint32_t oldr = vars.num_chunks;
    for (size_t i = 0; i < vars.D.size(); ++i) {
      const uint32_t r = vars.D[i].r;
      const uint32_t c = vars.D[i].c;
      tl->submit(new SparseProdTask(name, r, c, T(r), vars.D[i], H(c), r != oldr));
      oldr = r;
    }

    // F(r) = rhs(H(r), T(r))
    for (uint32_t r = 0; r < vars.num_chunks; ++r)
      tl->submit(new CalcRHSTask(name, vars.atm[r], H(r), T(r), vars.gh0, F(r), prio_row[r]));
  }
};

#endif // RHS_HPP_INCLUDED
