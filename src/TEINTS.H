#ifndef _TEINTS_H
#define _TEINTS_H

#include <stdio.h>

#include "cgf.h"
#include "rep.h"

double* build_two_electron_integrals(const CGF *cgfs, unsigned nrorbs);

void print_teints(const double* teints, unsigned nrorbs);

#endif