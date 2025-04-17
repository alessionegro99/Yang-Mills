#include "../include/measure.h"

void staple(double complex **restrict lattice,
            long int const *const restrict nnp,
            long int const *const restrict nnm, long int r, int dir,
            long int stvol, double complex *ris) {
  //               ^ dir
  //         (4)   |   (1)
  //     +----<----+---->----+
  //     |         |         |
  //  (5)|         |         |
  //     V         ^         V (2)
  //     |         |         |
  //     +---->----+----<----+---> i
  //     r1  (6)   r   (3)
  //

  int i, k;
  long int r1;
  double complex aux;

  *ris = 0.0 + I * 0.0;

  for (k = 1; k < STDIM; k++) {
    i = (dir + k) % STDIM;

    // forward
    aux = lattice[nnp[dirgeo(r, dir, stvol)]][i];        // 1
    aux *= conj(lattice[nnp[dirgeo(r, i, stvol)]][dir]); // 2
    aux *= conj(lattice[r][i]);                          // 3
    *ris += aux;

    // backward
    r1 = nnm[dirgeo(r, i, stvol)];
    aux = conj(lattice[nnp[dirgeo(r1, dir, stvol)]][i]); // 4
    aux *= conj(lattice[r1][dir]);
    aux *= lattice[r1][i];
    *ris += aux;
  }
}

double complex multihit(double complex **restrict lattice,
                        long int const *const restrict nnp,
                        long int const *const restrict nnm, long int r, int dir,
                        double epsilon, double beta, int nhits,
                        long int stvol) {
  int i, count;
  double oldS, newS, deltaS, aux;
  double complex stap, link, change, ris = 0.0 + 0.0 * I;

  link = lattice[r][dir];

  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  count = 0;
  for (i = 0; i < nhits; i++) {
    // metropolis
    aux = epsilon * (2.0 * myrand() - 1.0);
    change =
        1.0 / sqrt(1.0 + aux * aux) +
        I * aux /
            sqrt(1.0 + aux * aux); // complex phase that becomes 1 if epsilon=0

    oldS = -beta * creal(link * stap);
    newS = -beta * creal(link * stap * change);
    deltaS = newS - oldS;

    if (deltaS < 0) {
      link *= change;
    } else {
      if (myrand() < exp(-deltaS)) {
        link *= change;
      }
    }
    ris += link;
    count++;

    // overrelaxation
    if (cabs(stap) > 1.0e-10) {
      change = conj(link * stap * stap);
      change /= (cabs(stap) * cabs(stap));

      link = change;

      ris += link;
      count++;
    }
  }

  ris /= (double complex)(count);

  return ris;
}

double Wilsonloop(double complex **restrict lattice,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm, int Wt, int Ws,
                  double epsilon, double beta, long int stvol) {
  int i, dir;
  const int nhits = 10;
  long int r;
  double ris = 0.0;
  double complex aux;

  for (r = 0; r < stvol; r++) {
    for (dir = 1; dir < STDIM; dir++) {
      //      ^ 0
      //      |
      //   r1 +--------+r2
      //      |        |
      //      |Wt      |
      //      |        |
      //      |        |
      //      +--------+--> dir
      //     r    Ws   r3

      aux = 1.0 + 0.0 * I;

      if (Wt == 1 || Ws == 1) {
        for (i = 0; i < Wt; i++) {
          aux *= lattice[r][0];
          r = nnp[dirgeo(r, 0, stvol)];
        }
        // now we are in r1

        for (i = 0; i < Ws; i++) {
          aux *= lattice[r][dir];
          r = nnp[dirgeo(r, dir, stvol)];
        }
        // now we are in r2

        for (i = 0; i < Wt; i++) {
          r = nnm[dirgeo(r, 0, stvol)];
          aux *= conj(lattice[r][0]);
        }
        // now we are in r3

        for (i = 0; i < Ws; i++) {
          r = nnm[dirgeo(r, dir, stvol)];
          aux *= conj(lattice[r][dir]);
        }
      } else {
        for (i = 0; i < Wt; i++) {
          if (i == 0) {
            aux *= lattice[r][0];
          } else {
            aux *=
                multihit(lattice, nnp, nnm, r, 0, epsilon, beta, nhits, stvol);
          }
          r = nnp[dirgeo(r, 0, stvol)];
        }
        // now we are in r1

        for (i = 0; i < Ws; i++) {
          if (i == 0) {
            aux *= lattice[r][dir];
          } else {
            aux *= multihit(lattice, nnp, nnm, r, dir, epsilon, beta, nhits,
                            stvol);
          }
          r = nnp[dirgeo(r, dir, stvol)];
        }
        // now we are in r2

        for (i = 0; i < Wt; i++) {
          r = nnm[dirgeo(r, 0, stvol)];
          if (i == 0) {
            aux *= conj(lattice[r][0]);
          } else {
            aux *= conj(
                multihit(lattice, nnp, nnm, r, 0, epsilon, beta, nhits, stvol));
          }
        }
        // now we are in r3

        for (i = 0; i < Ws; i++) {
          r = nnm[dirgeo(r, dir, stvol)];
          if (i == 0) {
            aux *= conj(lattice[r][dir]);
          } else {
            aux *= conj(multihit(lattice, nnp, nnm, r, dir, epsilon, beta,
                                 nhits, stvol));
          }
        }
      }
      ris += creal(aux);
    }
  }

  ris /= (double)stvol;
  ris /= (double)(STDIM - 1);

  return ris;
}

double nonplanarWilsonloop(double complex **restrict lattice,
                           long int const *const restrict nnp,
                           long int const *const restrict nnm, int Wt, int Ws1,
                           int Ws2, double epsilon, double beta, long int stvol,
                           int dir1, int dir2) {
  int i;
  const int nhits = 10;
  long int r;
  double ris = 0.0;
  double complex aux;

  for (r = 0; r < stvol; r++) {
    //      ^ 0
    //      |
    //   r1 +--------+r2
    //      |        |
    //      |Wt      |
    //      |        |
    //      |        |
    //      +--------+--> dir1
    //     r    Ws1  r5

    //      ^ 0
    //      |
    //   r2 +--------+r3
    //      |        |
    //      |Wt      |
    //      |        |
    //      |        |
    //      +--------+--> dir2
    //     r5   Ws2  r4

    aux = 1.0 + 0.0 * I;

    if (Wt == 1 || Ws1 == 1 || Ws2 == 1) {
      for (i = 0; i < Wt; i++) {
        aux *= lattice[r][0];
        r = nnp[dirgeo(r, 0, stvol)];
      }
      // now we are in r1

      for (i = 0; i < Ws1; i++) {
        aux *= lattice[r][dir1];
        r = nnp[dirgeo(r, dir1, stvol)];
      }
      // now we are in r2

      for (i = 0; i < Ws2; i++) {
        aux *= lattice[r][dir2];
        r = nnp[dirgeo(r, dir2, stvol)];
      }
      // now we are in r3

      for (i = 0; i < Wt; i++) {
        r = nnm[dirgeo(r, 0, stvol)];
        aux *= conj(lattice[r][0]);
      }
      // now we are in r4

      for (i = 0; i < Ws2; i++) {
        r = nnm[dirgeo(r, dir2, stvol)];
        aux *= conj(lattice[r][dir2]);
      }
      // now we are in r5

      for (i = 0; i < Ws1; i++) {
        r = nnm[dirgeo(r, dir1, stvol)];
        aux *= conj(lattice[r][dir1]);
      }
      // now we are in r
    } else {
      for (i = 0; i < Wt; i++) {
        if (i == 0) {
          aux *= lattice[r][0];
        } else {
          aux *= multihit(lattice, nnp, nnm, r, 0, epsilon, beta, nhits, stvol);
        }
        r = nnp[dirgeo(r, 0, stvol)];
      }
      // now we are in r1

      for (i = 0; i < Ws1; i++) {
        if (i == 0) {
          aux *= lattice[r][dir1];
        } else {
          aux *=
              multihit(lattice, nnp, nnm, r, dir1, epsilon, beta, nhits, stvol);
        }
        r = nnp[dirgeo(r, dir1, stvol)];
      }
      // now we are in r2

      for (i = 0; i < Ws2; i++) {
        if (i == 0) {
          aux *= lattice[r][dir2];
        } else {
          aux *=
              multihit(lattice, nnp, nnm, r, dir2, epsilon, beta, nhits, stvol);
        }
        r = nnp[dirgeo(r, dir2, stvol)];
      }
      // now we are in r3

      for (i = 0; i < Wt; i++) {
        r = nnm[dirgeo(r, 0, stvol)];
        if (i == 0) {
          aux *= conj(lattice[r][0]);
        } else {
          aux *= conj(
              multihit(lattice, nnp, nnm, r, 0, epsilon, beta, nhits, stvol));
        }
      }
      // now we are in r4

      for (i = 0; i < Ws2; i++) {
        r = nnm[dirgeo(r, dir2, stvol)];
        if (i == 0) {
          aux *= conj(lattice[r][dir2]);
        } else {
          aux *= conj(multihit(lattice, nnp, nnm, r, dir2, epsilon, beta, nhits,
                               stvol));
        }
      }
      // now we are in r5

      for (i = 0; i < Ws1; i++) {
        r = nnm[dirgeo(r, dir1, stvol)];
        if (i == 0) {
          aux *= conj(lattice[r][dir1]);
        } else {
          aux *= conj(multihit(lattice, nnp, nnm, r, dir1, epsilon, beta, nhits,
                               stvol));
        }
      }
    }
    ris += creal(aux);
  }

  ris /= (double)stvol;
  // ris /= (double)(STDIM - 1);

  return ris;
}