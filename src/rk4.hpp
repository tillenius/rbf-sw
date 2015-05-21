#ifndef RK4_HPP_INCLUDED
#define RK4_HPP_INCLUDED

#include "rhs.hpp"
#include "addtask.hpp"
#include "steptask.hpp"

//=============================================================================
// RK4
//=============================================================================
struct RK4 {
  RHS &rhs;
  const uint32_t chunk_size;
  const uint32_t num_chunks;
  const uint32_t nd;

  BlockedVectorHandle<vec4> K[3], d[4];

  RK4(RHS &rhs_)
  : rhs(rhs_),
    chunk_size(rhs.vars.chunk_size),
    num_chunks(rhs.vars.num_chunks), 
    nd(rhs.vars.nd)
  {
    for (size_t i = 0; i < 3; ++i)
      K[i].init(*tl, rhs.vars, "K");
    for (size_t i = 0; i < 4; ++i)
      d[i].init(*tl, rhs.vars, "d");
  }

  void add(double time,
           BlockedVectorHandle<vec4> &dst,
           BlockedVectorHandle<vec4> &A, const double s,
           BlockedVectorHandle<vec4> &B) {
    for (uint32_t r = 0; r < num_chunks; ++r)
      tl->submit(new AdditionTask(time, dst(r), A(r), s, B(r), r));
  }

  void step(double time, BlockedVectorHandle<vec4> &Hdst, BlockedVectorHandle<vec4> &Hsrc, const double s) {
    for (uint32_t r = 0; r < num_chunks; ++r)
      tl->submit(new StepTask(time, Hdst(r), Hsrc(r), s, d[0](r), d[1](r), d[2](r), d[3](r), r));
  }

  template<typename TaskGenerator>
  void operator()(double time, double dt, BlockedVectorHandle<vec4> &H,
                  TaskGenerator &task_gen) {

    rhs.calc_rhs(time, time, H, d[0]);           // d1 = dt*rhs(t, H);

    add(time, K[0], H, 0.5*dt, d[0]);            // K0 = H + 0.5*d1;
    rhs.calc_rhs(time, time+0.5*dt, K[0], d[1]); // d2 = dt*rhs(t+0.5*dt, K0);

    add(time, K[1], H, 0.5*dt, d[1]);            // K1 = H + 0.5*d2;
    rhs.calc_rhs(time, time+0.5*dt, K[1], d[2]); // d3 = dt*rhs(t+0.5*dt, K1);

    add(time, K[2], H, dt, d[2]);                // K2 = H + d3;
    rhs.calc_rhs(time, time+dt, K[2], d[3]);     // d4 = dt*rhs(t+dt, K2);

    // H = H + 1.0/6.0*(d1 + 2.0*d2 + 2.0*d3 + d4);
    step(time, H, H, dt/6.0);

    task_gen.callback(H);
  }
};

#endif // RK4_HPP_INCLUDED
