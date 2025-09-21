#ifndef _GTO_H
#define _GTO_H

#include <math.h>
#include "mymath.h"
#include "gamma.h"

/*
 * Gaussian Type Orbital
 *
 * x^l * y^m * z^n * exp(-alpha * r)
 */
struct GTO {
    int l;    /* power in x */
    int m;    /* power in y */
    int n;    /* power in z */
    double alpha;
    double px;
    double py;
    double pz;
};

/*
 * Calculate normalization constant of GTO
 */
double gto_norm(const struct GTO* g);

#endif