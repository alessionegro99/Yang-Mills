#ifndef MEASURE_H
#define MEASURE_H

#include "geometry.h"
#include "random.h"
#include "macro.h"

void staple(double complex **restrict lattice,
            long int const *const restrict nnp,
            long int const *const restrict nnm, long int r, int dir,
            long int stvol, double complex *ris);

double complex multihit(double complex **restrict lattice,
                        long int const *const restrict nnp,
                        long int const *const restrict nnm, long int r, int dir,
                        double epsilon, double beta, int nhits, long int stvol);

double Wilsonloop(double complex **restrict lattice,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm, int Wt, int Ws,
                  double epsilon, double beta, long int stvol);

double nonplanarWilsonloop(double complex **restrict lattice,
                           long int const *const restrict nnp,
                           long int const *const restrict nnm, int Wt, int Ws1,
                           int Ws2, double epsilon, double beta, long int stvol,
                           int dir1, int dir2);

#endif