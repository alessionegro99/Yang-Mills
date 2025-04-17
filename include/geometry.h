#ifndef GEOMETRY_H
#define GEOMETRY_H

// cartesian coordinates -> lexicographic index
// lattice Nt*Ns^{stdim-1}
void cart_to_lex_st(long int *restrict lex, int const *const restrict cartcoord,
                    int Nt, int Ns, int stdim);

// lexicographic index -> cartesian coordinates
// lattice Nt*Ns^{stdim-1}
void lex_to_cart_st(int *restrict cartcoord, long int lex, int Nt, int Ns,
                    int stdim);

// cartesian coordinates -> lexicographic even/odd index
void cart_to_lexeo_st(long int *restrict lexeo,
                      int const *const restrict cartcoord, int Nt, int Ns,
                      int stdim, long int volume);

// lexicographic eo index -> cartesian coordinates
void lexeo_to_cart(int *restrict cartcoord, long int lexeo, int Nt, int Ns,
                   int stdim, long int volume);

// nnp[dirgeo(r, i, volume)] is the neighbor of "r" in positive direction "i" on
// a lattice of volume "volume" nnm[dirgeo(r, i, volume)] is the neighbor of "r"
// in negative direction "i" on a lattice of volume "volume"
inline long int dirgeo(long int lex, int i, long int volume) {
  return i * volume + lex;
}

// initialize geometry
// nnp[dirgeo(r,i,volume)]= next neighbor in positive "i" direction of site r
// nnm[dirgeo(r,i,volume)]= next neighbor in negative "i" direction of site r
void init_neighbors_st(long int *restrict nnp, long int *restrict nnm, int Nt,
                       int Ns, int stdim);

void test_geometry(long int const *const restrict nnp,
                   long int const *const restrict nnm, int Nt, int Ns,
                   int stdim);

#endif
