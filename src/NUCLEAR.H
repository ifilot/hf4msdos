#ifndef _NUCLEAR_H
#define _NUCLEAR_H

#include "cgf.h"

/*
 * Calculate nuclear integral for two CGFs
 */
double cgf_nuclear(const CGF& cgf1, const CGF& cgf2,
                   double cx, double cy, double cz);

/*
 * Calculate nuclear integral for two GTOs
 */
double gto_nuclear(const struct GTO *g1, const struct GTO *g2,
                   double nx, double ny, double nz);

/*
 * Build A array for nuclear integral. Allocates memory to populate result
 * vector. Note that allocated memory needs to be free'd later.
 */
void build_A_array(double *arr, unsigned imax,
                   int l1, int l2,
                   double pa, double pb, double cp, double g);

/*
 * Calculates single term of A array
 */
double a_term(int i, int r, int u, int l1, int l2,
              double pax, double pbx, double cpx, double gamma);

#endif