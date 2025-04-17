#include "../include/geometry.h"
#include "../include/macro.h"
#include "../include/random.h"
#include "../include/update.h"

int main(int argc, char **argv) {
  int i, Nt, Ns, Wt, Ws;
  double beta, rand, W;
  double complex **lattice;
  long int *nnp, *nnm;
  long int iter, sample, r, stvolume, acc, count;

  char datafile[STRING_LENGTH];
  FILE *fp;

  const int overrelax = 5;
  const int measevery = 50;
  const int unitarizeevery = 10;
  const double epsilon = 1;

  const unsigned long int seed1 = (unsigned long int)time(NULL);
  const unsigned long int seed2 = seed1 + 127;

  if (argc != 6) {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s Nt Ns beta sample datafile\n\n", argv[0]);
    fprintf(stdout, "  Nt = temporal size of the lattice\n");
    fprintf(stdout, "  Ns = spatial size of the lattice (space-time dimension "
                    "defined by macro STDIM)\n");
    fprintf(stdout, "  beta = coupling\n");
    fprintf(stdout, "  sample = number of drawn to be extracted\n");
    fprintf(stdout,
            "  datafile = name of the file on which to write the data\n\n");
    fprintf(stdout, "Compiled for:\n");
    fprintf(stdout, "  dimensionality = %d\n\n", STDIM);
    fprintf(stdout, "Output:\n");
    fprintf(stdout, "  Wilson loops (Wt, Ws) using the following order\n");
    fprintf(stdout, "  for(Ws=1; Ws<=Ns/4; Ws++){\n");
    fprintf(stdout, "      for(Wt=1; Wt<=MIN(Nt/4,8); Wt++){\n");
    fprintf(stdout, "         Ws, Wt Wilson loop }}\n");

    return EXIT_SUCCESS;
  } else {
    // read input values
    Nt = atoi(argv[1]);
    Ns = atoi(argv[2]);
    beta = atof(argv[3]);
    sample = atol(argv[4]);

    if (strlen(argv[5]) >= STRING_LENGTH) {
      fprintf(stderr,
              "File name too long. Increse STRING_LENGTH or shorten the name "
              "(%s, %d)\n",
              __FILE__, __LINE__);
      return EXIT_FAILURE;
    } else {
      strcpy(datafile, argv[5]);
    }
  }

  // initialize random number generator
  myrand_init(seed1, seed2);

  // compute the spacetime volume
  stvolume = Nt;
  for (i = 1; i < STDIM; i++) {
    stvolume *= Ns;
  }

  // allocate the lattice
  // and next neighbors: nnp[dirgeo(r, i, volume)]= next neighbor in positive
  // "i" direction of site r
  lattice = (double complex **)malloc((unsigned long int)(stvolume) *
                                      sizeof(double complex *));
  if (lattice == NULL) {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
  }
  for (r = 0; r < stvolume; r++) {
    lattice[r] = (double complex *)malloc((unsigned long int)(STDIM) *
                                          sizeof(double complex));
    if (lattice[r] == NULL) {
      fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
      return EXIT_FAILURE;
    }
  }
  nnp = (long int *)malloc((unsigned long int)(STDIM * stvolume) *
                           sizeof(long int));
  if (nnp == NULL) {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
  }
  nnm = (long int *)malloc((unsigned long int)(STDIM * stvolume) *
                           sizeof(long int));
  if (nnm == NULL) {
    fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  // initialize nnp and nnm
  init_neighbors_st(nnp, nnm, Nt, Ns, STDIM);

  // initialize lattice to ordered start
  for (r = 0; r < stvolume; r++) {
    for (i = 0; i < STDIM; i++) {
      lattice[r][i] = 1.0 + I * 0.0;
    }
  }

  // open data file
  fp = fopen(datafile, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile,
            __FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  acc = 0.0;
  count = 0;
  for (iter = 0; iter < sample; iter++) {
    rand = myrand();

    if (rand < 1.0 / (double)overrelax) {
      count++;
      // metropolis
      for (r = 0; r < stvolume; r++) {
        for (i = 0; i < STDIM; i++) {
          acc += metropolis(lattice, nnp, nnm, r, i, epsilon, beta, stvolume);
        }
      }
    } else {
      // overrelaxation
      for (r = 0; r < stvolume; r++) {
        for (i = 0; i < STDIM; i++) {
          overrelaxation(lattice, nnp, nnm, r, i, stvolume);
        }
      }
    }

    // reunitarize
    if (iter % unitarizeevery == 0) {
      for (r = 0; r < stvolume; r++) {
        for (i = 0; i < STDIM; i++) {
          lattice[r][i] /= cabs(lattice[r][i]);
        }
      }
    }

    // perform measures
    if (iter % measevery == 0) {
      // perform measures
      for (Ws = 1; Ws <= Ns / 4; Ws++) {
        for (Wt = 1; Wt <= MIN(Nt / 4, 8); Wt++) {
          W = Wilsonloop(lattice, nnp, nnm, Wt, Ws, epsilon, beta, stvolume);
          fprintf(fp, "%.12f ", W);
        }
      }
      fprintf(fp, "\n");
    }

    // perform measures
    if (iter % measevery == 0) {
      // perform measures

      fprintf(fp, "%.12f ",
              (nonplanarWilsonloop(lattice, nnp, nnm, 1, 1, 1, epsilon, beta,
                                   stvolume, 1, 2) +
               nonplanarWilsonloop(lattice, nnp, nnm, 1, 1, 1, epsilon, beta,
                                   stvolume, 2, 1)) /
                  2.0);

      fprintf(fp, "\n");
    }
  }

  printf("Acceptance rate %f\n",
         (double)acc / (double)count / (double)stvolume / (double)STDIM);

  // close datafile
  fclose(fp);

  // free memory
  for (r = 0; r < stvolume; r++) {
    free(lattice[r]);
  }
  free(lattice);
  free(nnp);
  free(nnm);

  return EXIT_SUCCESS;
}
