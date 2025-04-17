#ifndef UPDATE_H
#define UPDATE_H

#include "measure.h"

int metropolis(double complex **restrict lattice,
               long int const *const restrict nnp,
               long int const *const restrict nnm, long int r, int dir,
               double epsilon, double beta, long int stvol);

void overrelaxation(double complex **restrict lattice,
                    long int const *const restrict nnp,
                    long int const *const restrict nnm, long int r, int dir,
                    long int stvol);

#endif