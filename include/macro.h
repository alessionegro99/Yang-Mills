#ifndef MACRO_H
#define MACRO_H

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// when using C99 M_PI is not defined in math.h header!
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b)) // min of two numbers

#define STDIM 3 // space-time dimensionality
#define STRING_LENGTH 50

#endif