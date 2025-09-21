#include "gto.h"
#include "integrals.h"

/*
 * Calculate normalization constant of GTO
 */
double gto_norm(const struct GTO *g) {
    double ps = (double)(g->l + g->m + g->n);
    double nom = pow(2.0, 2.0 * ps + 3.0 / 2.0) *
                 pow(g->alpha, ps + 3.0 / 2.0);
    double denom = fact2(2 * g->l-1) *
                   fact2(2 * g->m-1) *
                   fact2(2 * g->n-1) *
                   pow(M_PI, 3.0 / 2.0);
    return sqrt(nom / denom);
}
