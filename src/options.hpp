#ifndef OPTIONS_HPP_INCLUDED
#define OPTIONS_HPP_INCLUDED

#include "sg/option/instr_trace.hpp"
#include "sg/option/taskqueue_priopinned.hpp"

#ifdef USE_MPI

#ifdef DEBUG_TRACE
#include "sgmpi/option/mpidebugdefaults.hpp"
struct Options : public sgmpi::DebugOptions<Options> {
#else
struct Options : public sgmpi::DefaultOptions<Options> {
#endif

#else
struct Options : public DefaultOptions<Options> {
#endif

#ifdef DEBUG_TRACE
  typedef Trace<Options> Instrumentation;
#endif

  typedef Enable TaskName;
  typedef Enable PassTaskExecutor; // Enable for non-MPI SuperGlue to have common interface
  typedef TaskQueuePrioPinned<Options> WaitListType;
  typedef TaskQueuePrioPinned<Options> ReadyListType;
};

#endif // OPTIONS_HPP_INCLUDED
