#ifndef OPTIONS_HPP_INCLUDED
#define OPTIONS_HPP_INCLUDED

#include "sg/option/taskqueue_priopinned.hpp"
#include "sg/option/log2.hpp"

template<typename Options>
struct MyTaskBase : public TaskBaseDefault<Options> {
  double time;
  MyTaskBase() : time(-1) {}
};

template<typename Options>
struct MyTrace {
  Time::TimeUnit start, stop;
  MyTrace(int threadid) {
    Log2<Options, double>::register_thread(threadid);
  }
  void run_task_before(sg::TaskBase<Options> *) {
    start = Time::getTime();
  }
  void run_task_after(sg::TaskBase<Options> *task) {
    stop = Time::getTime();
    Log2<Options, double>::log(task->time, start, stop);
  }
  static void dump(const char *name) {
    Log2<Options, double>::dump(name);
  }
};

template<typename Options>
class MyMPITrace {
public:
  struct MPIInstrData {
    Time::TimeUnit mpiinstr_start;
  };

  static void start() {}
  template<typename Req>
  static void start(Req &req) {
    req.mpiinstr_start = Time::getTime();
  }
  template<typename Req>
  static void recv(Req &req, Time::TimeUnit stop) {
    Log2<Options, double>::log(-2.0, stop, stop);
  }
  template<typename Req>
  static void send(Req &req, Time::TimeUnit stop) {
    Log2<Options, double>::log(-1.0, req.mpiinstr_start, stop);
  }
};

#ifdef USE_MPI
  #define RBFSW_OPTIONS_BASE sgmpi::DefaultOptions<Options>
#else
  #define RBFSW_OPTIONS_BASE sg::DefaultOptions<Options>
#endif

struct Options : public RBFSW_OPTIONS_BASE {
#ifdef DEBUG_TRACE
  typedef MyTrace<Options> Instrumentation;
  #ifdef USE_MPI
    typedef MyMPITrace<Options> MPIInstrumentation;
  #endif
#endif
  typedef MyTaskBase<Options> TaskBaseType;

  typedef Enable PassTaskExecutor; // Enable for non-MPI SuperGlue to have common interface
  typedef TaskQueuePrioPinned<Options> WaitListType;
  typedef TaskQueuePrioPinned<Options> ReadyListType;
};

#undef RBFSW_OPTIONS_BASE

#endif // OPTIONS_HPP_INCLUDED
