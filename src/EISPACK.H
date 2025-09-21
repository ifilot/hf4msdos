#ifndef _EISPACK_H
#define _EISPACK_H

#include "math.h"

double r8_abs(double x);
double r8_max(double x, double y);
double r8_min(double x, double y);
double r8_sign(double x);
double r8_epsilon();
double pythag(double a, double b);

void tred1(int n, double a[], double d[], double e[], double e2[]);
void tred2(int n, double a[], double d[], double e[], double z[]);
int  tql2(int n, double d[], double e[], double z[]);
int  tqlrat(int n, double d[], double e2[]);

int  rs(int n, double a[], double w[], int matz, double z[]);

#endif