#ifndef _OVERLAP_H
#define _OVERLAP_H

#include "cgf.h"

/*
 * Calculate the overlap integral between two CGF objects
 */
double cgf_overlap(const CGF& cgf1, const CGF& cgf2);

/*
 * Calculate overlap integral of two GTOs
 */
double gto_overlap(const struct GTO *g1, const struct GTO *g2);

/*
 * Calculate overlap integral of two GTOs with modified Cartesian powers.
 * This function is used as a subroutine in the kinetic overlap integral
 */
double gto_overlap_modcoeff(const struct GTO *g1, const struct GTO *g2,
			    int l1, int m1, int n1,
			    int l2, int m2, int n2);

/*
 * Calculate 1D overlap integral for a GTO
 */
double gto_overlap_1d(int l1, int l2, double x1, double x2, double gamma);

#endif