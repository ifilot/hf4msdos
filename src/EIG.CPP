#include "eig.h"

void eigh(const Matrix& A, Matrix *V, Matrix *D) {
    unsigned int n = A.cols();     // size of the matrix
    double *w = new double[n];     // store eigenvalues
    int i,j;

    static int matz = 1;           // return both eigenvalues as well as eigenvectors

    unsigned int ierr = 0;
    Matrix Acopy = A;              // changes can be made to A, so create copy
    ierr = rs(n, Acopy.data(), w, matz, V->data());

    for(i=0; i<n; i++) {
        D->at(i,i) = w[i];
    }

    V->self_transpose();

    delete[] w;
}