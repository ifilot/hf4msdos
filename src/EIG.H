#ifndef _EIG_H
#define _EIG_H

#include <math.h>
#include <conio.h>

#include "matrix.h"
#include "eispack.h"

void eigh(const Matrix& A, Matrix *V, Matrix *D);

#endif