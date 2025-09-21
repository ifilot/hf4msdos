#ifndef _REPULSION_H
#define _REPULSION_H

#include "cgf.h"
#include "gamma.h"
#include "mymath.h"

double cgf_repulsion(const CGF& cgf1, const CGF& cgf2,
                     const CGF& cgf3, const CGF& cgf4);

double gto_repulsion(const struct GTO *g1, const struct GTO *g2, 
                     const struct GTO *g3, const struct GTO *g4);

void build_B_array(double *arr, int l1,int l2,int l3,int l4,
                   double p, double a, double b, double q, double c, double d,
                   double g1, double g2, double delta);

double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2, int l3, int l4,
double px, double ax, double bx, double qx, double cx, double dx, double gamma1,
double gamma2, double delta);

double fB(int i, int l1, int l2, double p, double a, double b, int r, double q);
double B0(int i, int r, double q);
double fact_ratio2(int a, int b);

const unsigned teindex(unsigned i, unsigned j, unsigned k, unsigned l);

#endif //_REPULSION_H