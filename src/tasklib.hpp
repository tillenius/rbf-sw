#ifndef TASKLIB_HPP_INCLUDED
#define TASKLIB_HPP_INCLUDED

#ifdef USE_MPI
extern sgmpi::MPISuperGlue<Options> *tl;
#else
extern SuperGlue<Options> *tl;

namespace sgmpi {

template<typename Options>
struct MPIHandle : public Handle<Options> {};

template<typename Options>
struct MPITask : public Task<Options> {};

} // namespace sgmpi

#endif

#endif // TASKLIB_HPP_INCLUDED
