#ifdef USE_MPI
#include "sgmpi/superglue_mpi.hpp"
#include "options.hpp"
typedef sgmpi::MPISuperGlue<Options> tasklib_t;
#else
#include "sg/superglue.hpp"
#include "options.hpp"
typedef SuperGlue<Options> tasklib_t;
#endif

#include "matrix.hpp"
#include "matrixutil.hpp"
#include "blocked_vec.hpp"
#include "variables.hpp"
#include "rhs.hpp"
#include "rk4.hpp"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <vector>

typedef unsigned int uint32_t;

tasklib_t *tl;
variables vars;
std::vector<uint32_t> prio_row;

std::string g_mldatapath;
uint32_t g_chunksize;
double g_dt, g_endtime;
double g_continuetime = -1.0;

//=============================================================================
// verify
//=============================================================================
std::string verify(BlockedVectorHandle<vec4> &H) {
  std::stringstream ss;
  double u,v,w,h,umin,vmin,wmin,hmin,umax,vmax,wmax,hmax;
  umin = umax = u = *H(0)(0).elem(0);
  vmin = vmax = v = *H(0)(0).elem(1);
  wmin = wmax = w = *H(0)(0).elem(2);
  hmin = hmax = h = *H(0)(0).elem(3);

  for (size_t i = 0; i < H.size(); ++i) {
    for (size_t j = 0; j < H(i).size(); ++j) {
      double tu = *H(i)(j).elem(0); u += tu; umin = std::min(umin, u); umax = std::max(umax, u);
      double tv = *H(i)(j).elem(1); v += tv; vmin = std::min(vmin, v); vmax = std::max(vmax, v);
      double tw = *H(i)(j).elem(2); w += tw; wmin = std::min(wmin, w); wmax = std::max(wmax, w);
      double th = *H(i)(j).elem(3); h += th; hmin = std::min(hmin, h); hmax = std::max(hmax, h);
    }
  }
  ss << " dbg=[ " << u << " " << v << " " << w << " " << h << " "
     << umin << " " << vmin << " " << wmin << " " << hmin << " "
     << umax << " " << vmax << " " << wmax << " " << hmax << " ]";
  return ss.str();
}

struct TaskGenerator;

struct SaveTask : public sgmpi::MPITask<Options> {
  std::string filename;
  VectorBlockHandle<vec4> &H;
  bool overwrite;
  SaveTask(const char *filename_, VectorBlockHandle<vec4> &H_, sgmpi::MPIHandle<Options> &filehandle, bool o_)
  : filename(filename_), H(H_), overwrite(o_) {
    this->time = -4.0;
    register_access(ReadWriteAdd::read, H.handle);
    register_access(ReadWriteAdd::write, filehandle);
  }
  void run(TaskExecutor<Options> &) {
    FILE *fout = fopen(filename.c_str(), overwrite ? "w" : "a");
    if (fout == NULL) { std::cerr<<"Error opening file " << filename << std::endl; exit(1); }
    if (overwrite) {
      size_t rows = 4;
      fwrite(&rows, 1, sizeof(size_t), fout);
      size_t cols = vars.nd;
      fwrite(&cols, 1, sizeof(size_t), fout);
    }
    fwrite(H.get_data(), H.size(), sizeof(vec4), fout);
    fclose(fout);
  }
};

struct GenTaskCallback : public Task<Options> {
  TaskGenerator &tg;
  const VectorBlockHandle<vec4> *block;
  GenTaskCallback(TaskGenerator &tg_, BlockedVectorHandle<vec4> &H);
  void run(TaskExecutor<Options> &);
};

#ifdef RBFSW_DEBUG
const double save_interval = 60.0*60.0*3.0; // save every 3 hours
#endif

struct TaskGenerator {
  double t;
  double dt;
  double t_end;
#ifdef RBFSW_DEBUG
  double next_save;
  SpinLock abort_spinlock;
  volatile bool aborted;
#endif
  Handle<Options> handle; // to disallow several concurrent TaskGenerator

  RHS rhs;
  RK4 rk4;

  BlockedVectorHandle<vec4> &H;

  sgmpi::MPIHandle<Options> filehandle;

  TaskGenerator(variables &vars,
                double dt_, double t_end_, BlockedVectorHandle<vec4> &H_)
  : t(g_continuetime != -1.0 ? g_continuetime : 0.0),
    dt(dt_), t_end(t_end_),
#ifdef RBFSW_DEBUG
    next_save(t+save_interval), 
    aborted(false),
#endif
    rhs(vars), rk4(rhs), H(H_)
  {}

  void callback(BlockedVectorHandle<vec4> &H);

#ifdef RBFSW_DEBUG
  void generate(bool abort_now) {
#else
  void generate() {
#endif
    if (t >= t_end) {
      //fprintf(stderr, "Generate %f %f -- nope\n", t, t_end);
      return;
    }

#ifdef RBFSW_DEBUG
    if (aborted)
      return;

    if (abort_now) {
      SpinLockScoped hold(abort_spinlock);
      if (aborted)
        return;
      aborted = true;
      printf("Abort %f / %f\n", t, t_end);
      return;
    }

    if (t >= next_save) {
      next_save += save_interval;
      char fname[80];
      sprintf(fname, "results-t%d.txt", (int) t);
      for (size_t i = 0; i < H.size(); ++i)
        tl->submit(new SaveTask(fname, H(i), filehandle, i == 0));
    }
#endif

    rk4(t, dt, H, *this);
    t += dt;
  }
};

GenTaskCallback::GenTaskCallback(TaskGenerator &tg_, BlockedVectorHandle<vec4> &H)
: tg(tg_) {
  this->time = -5.0;
  register_access(ReadWriteAdd::write, tg_.handle);
#ifdef USE_MPI
  block = &H.register_access_any(tl->get_rank(), *this, ReadWriteAdd::read, 0); // TODO: don't require all tasks to finish!
#else
  register_access(ReadWriteAdd::read, H[0].handle);
  block = &H[0];
#endif
}

void GenTaskCallback::run(TaskExecutor<Options> &) {
#ifdef RBFSW_DEBUG
  tg.generate( isnan(*(double *) block->get_data()) );
#else
  tg.generate();
#endif
}

void TaskGenerator::callback(BlockedVectorHandle<vec4> &H) {
  tl->submit(new GenTaskCallback(*this, H));
}
//=============================================================================
// shallow water simulation
//=============================================================================
void sw(tasklib_t &tl_) {

  tl = &tl_;

#ifdef USE_MPI
  int rank = tl->get_rank();
  int num_ranks = tl->get_num_ranks();
  int num_cores = tl->get_num_cpus();
#else
  int rank = 0;
  int num_ranks = 1;
  int num_cores = tl->get_num_cpus();
#endif

//  fprintf(stderr, "# Load From File \n");
  vars.init(*tl, g_mldatapath.c_str(), num_cores, num_ranks, rank, g_chunksize);

  BlockedVectorHandle<vec4> H;
  if (g_continuetime == -1.0)
    vars.load(*tl, (g_mldatapath + "/H").c_str(), H);
  else {
    char fname[80];
    sprintf(fname, "results-t%d.txt", (int) g_continuetime);
    vars.load(*tl, fname, H);
    if (rank == 0)
      std::cout << "... continue from " << g_continuetime << std::endl;
  }

  if (rank == 0) {
    std::cout
       << "chunk_size= " << vars.chunk_size
       << " num_chunks= " << vars.num_chunks
       << " dt= " << g_dt
       << " end= " << g_endtime
       << " data= " << g_mldatapath
       << " cpus= " << num_cores << std::endl;
  }

  // mark blocks in H that will be communicated
#ifdef USE_MPI
  prio_row = std::vector<uint32_t>(vars.num_chunks, false);
  for (size_t i = 0; i < vars.D.size(); ++i) {
    const uint32_t r = vars.D[i].r;
    const uint32_t c = vars.D[i].c;
    if (H(c).handle.get_rank() != H(r).handle.get_rank())
      prio_row[r] = true;
  }
#endif // USE_MPI

  TaskGenerator task_gen(vars, g_dt, g_endtime, H);

#ifdef USE_MPI
  tl->mpi_barrier();
#endif // USE_MPI

  Time::TimeUnit start = Time::getTime();

  for (size_t i = 0; i < 5; ++i)
    task_gen.callback(H);

  tl->barrier();

  Time::TimeUnit stop = Time::getTime();

  tl->wait(task_gen.handle);

//  fprintf(stderr, "### DONE (%llu)\n", stop-start);

  if (rank == 0) {
    std::stringstream ss;
    ss << "size = " << H.size() 
       << " chunk_size= " << vars.chunk_size
       << " num_chunks= " << vars.num_chunks
       << " cpus= " << num_cores 
       << " time= " << std::setfill('0') << std::setw(1) << (stop - start) / 1000000 << "." 
       << std::setfill('0') << std::setw(6) << (stop - start) % ((stop - start) / 1000000);

    ss << verify(H);

    ss << std::endl;
    printf("%s", ss.str().c_str());

    const char *filename = "result.txt";
    FILE *fout = fopen(filename, "w");
    if (fout == NULL) { std::cerr<<"Error opening file " << filename << std::endl; exit(1); }
    size_t rows = 4;
    fwrite(&rows, 1, sizeof(size_t), fout);
    size_t cols = vars.nd;
    fwrite(&cols, 1, sizeof(size_t), fout);
    for (uint32_t r = 0; r < H.size(); ++r)
      fwrite(H(r).get_data(), H(r).size(), sizeof(vec4), fout);
    fclose(fout);
  }

  {
    std::stringstream ss;
    ss << "trace-" << rank << ".trace";
    Log2<Options, double>::dump(ss.str().c_str(), rank);
  }
}

//=============================================================================
// main
//=============================================================================
int main(int argc, char *argv[]) {

  assert(sizeof(uint32_t) == 4);

  if (argc != 5 && argc != 6) {
    fprintf(stderr, "Usage: %s <mldatapath> <chunk_size> <time_step> <end_time> [continue_from]\n", argv[0]);
    return 1;
  }

  g_mldatapath = argv[1];
  g_chunksize = atoi(argv[2]);
  g_dt = atof(argv[3]);
  g_endtime = atof(argv[4]);
  if (argc == 6)
    g_continuetime = atof(argv[5]);

#ifdef USE_MPI
  sgmpi::mpisuperglue<Options>(sw);
#else
  SuperGlue<Options> tm;
  sw(tm);
#endif
  return 0;
}
