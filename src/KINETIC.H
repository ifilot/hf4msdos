#ifndef _KINETIC_H
#define _KINETIC_H

#include "overlap.h"

/*
 * Calculate the kinetic integral between two CGF objects
 */
double cgf_kinetic(const CGF& cgf1, const CGF& cgf2);

/*
 * Calculate kinetic integral of two GTOs
 */
double gto_kinetic(const struct GTO *g1, const struct GTO *g2);

#endif