/*
 * References:
 * [1] Numerical Recipes
 * William H. Press, Saul A. Teukolsky, William T., 
 * Vetterling and Brian P. Flannery
 * 3rd edition page 261
 * ISBN-13: 978-0521880688
 */

#ifndef _GAMMA_H
#define _GAMMA_H

#include <iostream.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "mymath.h"

double fgamma(double m, double x);
double gamm_inc(double a, double x);
double gammp(double m, double x);
double gser(double a, double x);
double gammln(double xx);
double gcf(double a, double x);
double gammpapprox(double a, double x, int psig);

#endif