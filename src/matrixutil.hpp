#ifndef MATRIXUTIL_HPP_INCLUDED
#define MATRIXUTIL_HPP_INCLUDED

#include "quad.hpp"

template<typename value_t> void myfprintf(FILE *, value_t &);

template<>
void myfprintf(FILE *fout, vec4 &q) {
  const double *d = (double *) &q;
  fprintf(fout, "%f %f %f %f\n", d[0], d[1], d[2], d[3]);
}

template<>
void myfprintf(FILE *fout, tri<vec4> &q) {
  const double *d = (double *) &q;
  fprintf(fout, "%f %f %f %f %f %f %f %f %f %f %f %f\n",
    d[0], d[1], d[2], d[3], 
    d[4], d[5], d[6], d[7],
    d[8], d[9], d[10], d[11]);
}

template<typename value_t>
void saveText(const char *name, std::vector<value_t, AlignmentAllocator<value_t, 32> > data) {
  FILE *fout = fopen(name, "w");
  if (fout == NULL) { std::cerr<<"Error opening file " << name << std::endl; exit(1); }
  for (uint32_t r = 0; r < data.size(); ++r)
    myfprintf(fout, data[r]);
  fclose(fout);
}

#endif // MATRIXUTIL_HPP_INCLUDED
