#include "teints.h"

// build two-electron matrices
double* build_two_electron_integrals(const CGF *cgfs, unsigned nrorbs) {
    int i = 0; // top-level scope counter
    size_t nrtei = teindex(nrorbs, nrorbs, nrorbs, nrorbs) + 1;
    double* teints = (double*)malloc(nrtei * sizeof(double));

    for(i=0; i<nrtei; ++i) {
        teints[i] = -1.0;       // set some uninitialized value
    }

    // construct joblist
    int *jobcoeff = new int[4*nrtei];
    int *jobindices = new int[nrtei];

    // construct joblist
    size_t ctr = 0;
    for(i=0; i<nrorbs; ++i) {
        for(int j=0; j<=i; ++j) {
            int ij = i*(i+1)/2 + j;
            for(int k=0; k<nrorbs; ++k) {
                for(int l=0; l<=k; ++l) {
                    int kl = k*(k+1)/2 + l;
                    if(ij <= kl) {
                        int index = teindex(i,j,k,l);
                        if(teints[index] < 0) {
                            jobindices[ctr] = index;
                            int *ptr = &jobcoeff[ctr << 2];
                            *(ptr++) = i;
                            *(ptr++) = j;
                            *(ptr++) = k;
                            *(ptr++) = l;
                            ctr++;
                        }
                    }
                }
            }
        }
    }

    // build the integrals
    nrtei = ctr;
    printf("Evaluating %i two-electron integrals...\n", nrtei);
    for(i=0; i<ctr; ++i) {
        int index = jobindices[i];
        int *ptr = &jobcoeff[i << 2];
        int idx0 = *(ptr++);
        int idx1 = *(ptr++);
        int idx2 = *(ptr++);
        int idx3 = *ptr;
        teints[index] = cgf_repulsion(cgfs[idx0], cgfs[idx1], cgfs[idx2], cgfs[idx3]);
        printf("\r%i / %i integrals evaluated", i+1, nrtei);
    }
    cout << endl;

    // clean up the mess
    delete[] jobcoeff;
    delete[] jobindices;

    return teints;
}

void print_teints(const double* teints, unsigned nrorbs) {
    int i = 0; // top-level scope counter
    unsigned ctr = 0;

    // build the interals
    for(i=0; i<nrorbs; ++i) {
        for(int j=0; j<=i; ++j) {
            int ij = i*(i+1)/2 + j;
            for(int k=0; k<nrorbs; ++k) {
                for(int l=0; l<=k; ++l) {
                    int kl = k*(k+1)/2 + l;
                    if(ij <= kl) {
                        int index = teindex(i,j,k,l);
                        printf("(%i,%i,%i,%i) = %6.4f\n",
                                i+1,j+1,k+1,l+1,teints[index]);
                        ctr++;
                        if(ctr % 40 == 0) {
                            getch();
                        }
                    }
                }
            }
        }
    }
}