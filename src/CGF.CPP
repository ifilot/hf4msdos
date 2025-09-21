#include "cgf.h"

CGF::CGF() {}

/*
 * Create a new CGF
 *
 * @param N : number of GTOs (fixed)
 */
CGF::CGF(unsigned _N) {
    this->set_ngtos(_N);
}

void CGF::set_ngtos(unsigned _N) {
    this->N = _N;
    this->gto_list = (GTO*)malloc(this->N * sizeof(GTO));
    this->gto_coeff = (double*)calloc(this->N, sizeof(double));
}

/*
 * Set a GTO in the CGF
 *
 * @param idx : GTO index
 * @param coeff : GTO coefficient in CGF
 * @param l : power in x
 * @param m : power in y
 * @param n : power in z
 * @param alpha : exponent
 * @param px : x-position of GTO
 * @param py : y-position of GTO
 * @param pz : z-position of GTO
 */
void CGF::set_gto(unsigned idx,
                  double coeff,
                  int l,
                  int m,
                  int n,
                  double alpha,
                  double px,
                  double py,
                  double pz) {
    if(idx >= this->N) {
        cerr << "Invalid index!" << endl;
        exit(-1);
    }
    struct GTO gto = {l, m, n, alpha, px, py, pz};
    this->gto_list[idx] = gto;
    this->gto_coeff[idx] = gto_norm(&this->gto_list[idx]) * coeff;
}

void CGF::print() const {
    printf("CGF:\n");
    for(int i=0; i<this->N; ++i) {
        GTO *g = &this->gto_list[i];
        printf("  %i (%6.4f, %6.4f) (%i,%i,%i) (%6.4f,%6.4f,%6.4f)\n",
            i+1, this->gto_coeff[i], g->alpha,
            g->l, g->m, g->n,
            g->px, g->py, g->pz);
    }
}