#include "../include/update.h"

// return 1 if the update is accepted, 0 otherwise
int metropolis(double complex **restrict lattice,
               long int const *const restrict nnp,
               long int const *const restrict nnm, long int r, int dir,
               double epsilon, double beta, long int stvol) {
  int acc = 0;
  double oldS, newS, deltaS, aux;
  double complex stap, change;

  aux = epsilon * (2.0 * myrand() - 1.0);
  change =
      1.0 / sqrt(1.0 + aux * aux) +
      I * aux /
          sqrt(1.0 + aux * aux); // complex phase that becomes 1 if epsilon=0

  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  oldS = -beta * creal(lattice[r][dir] * stap);
  newS = -beta * creal(lattice[r][dir] * stap * change);
  deltaS = newS - oldS;

  if (deltaS < 0) {
    lattice[r][dir] *= change;
    acc = 1;
  } else {
    if (myrand() < exp(-deltaS)) {
      lattice[r][dir] *= change;
      acc = 1;
    }
  }

  return acc;
}

void overrelaxation(double complex **restrict lattice,
                    long int const *const restrict nnp,
                    long int const *const restrict nnm, long int r, int dir,
                    long int stvol) {
  double complex stap, aux;

  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  if (cabs(stap) > 1.0e-10) {
    aux = conj(lattice[r][dir] * stap * stap);
    aux /= (cabs(stap) * cabs(stap));

    lattice[r][dir] = aux;
  }
}

