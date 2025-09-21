#include "rep.h"

/*
 * Calculate the nuclear integral between two CGF objects
 */
double cgf_repulsion(const CGF& cgf1, const CGF& cgf2,
                     const CGF& cgf3, const CGF& cgf4) {
    double sum = 0;
    const unsigned char N1 = cgf1.get_number_gtos();
    const unsigned char N2 = cgf2.get_number_gtos();
    const unsigned char N3 = cgf1.get_number_gtos();
    const unsigned char N4 = cgf2.get_number_gtos();

    for(unsigned char i=0; i<N1; ++i) {
    	const struct GTO *g1 = cgf1.get_gto_ptr(i);
	    const double c1 = cgf1.get_coeff(i);
     	for(unsigned char j=0; j<N2; ++j) {
	        const struct GTO *g2 = cgf2.get_gto_ptr(j);
    	    const double c2 = cgf2.get_coeff(j);
            for(unsigned char j=0; j<N3; ++j) {
                const struct GTO *g3 = cgf3.get_gto_ptr(j);
	            const double c3 = cgf3.get_coeff(j);
                for(unsigned char j=0; j<N4; ++j) {
	                const struct GTO *g4 = cgf4.get_gto_ptr(j);
                    const double c4 = cgf4.get_coeff(j);
                    double term = c1 * c2 * c3 * c4 * gto_repulsion(g1, g2, g3, g4);
                    sum += term;
                }
            }
        }
    }

    return sum;
}

double gto_repulsion(const struct GTO *g1, const struct GTO *g2, 
                     const struct GTO *g3, const struct GTO *g4) {

    double dx = g1->px - g2->px;
    double dy = g1->py - g2->py;
    double dz = g1->pz - g2->pz;
    const double rab2 = (dx * dx + dy * dy + dz * dz);

    dx = g3->px - g4->px;
    dy = g3->py - g4->py;
    dz = g3->pz - g4->pz;
    const double rcd2 = (dx * dx + dy * dy + dz * dz);

    // calculate Gaussian product centers
    double gamma1 = g1->alpha + g2->alpha;
    const double px = (g1->alpha * g1->px + g2->alpha * g2->px) / gamma1;
    const double py = (g1->alpha * g1->py + g2->alpha * g2->py) / gamma1;
    const double pz = (g1->alpha * g1->pz + g2->alpha * g2->pz) / gamma1;

    double gamma2 = g3->alpha + g4->alpha;
    const double qx = (g3->alpha * g3->px + g4->alpha * g4->px) / gamma2;
    const double qy = (g3->alpha * g3->py + g4->alpha * g4->py) / gamma2;
    const double qz = (g3->alpha * g3->pz + g4->alpha * g4->pz) / gamma2;

    dx = px - qx;
    dy = py - qy;
    dz = pz - qz;
    double rpq2 = (dx * dx + dy * dy + dz * dz);

    double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);

    // build arrays and initialize to zero
    unsigned suml = g1->l + g2->l + g3->l + g4->l;
    unsigned summ = g1->m + g2->m + g3->m + g4->m;
    unsigned sumn = g1->n + g2->n + g3->n + g4->n;
    double *bx = new double[suml+1];
    memset(bx, 0.0, sizeof(double) * (suml+1));
    double *by = new double[summ+1];
    memset(by, 0.0, sizeof(double) * (summ+1));
    double *bz = new double[sumn+1];
    memset(bz, 0.0, sizeof(double) * (sumn+1));

    // populate arrays
    build_B_array(bx, g1->l, g2->l, g3->l, g4->l,
            px, g1->px, g2->px, qx, g3->px, g4->px,
            gamma1, gamma2, delta);
    build_B_array(by, g1->m, g2->m, g3->m, g4->m,
            py, g1->py, g2->py, qy, g3->py, g4->py,
            gamma1, gamma2, delta);
    build_B_array(bz, g1->n, g2->n, g3->n, g4->n,
            pz, g1->pz, g2->pz, qz, g3->pz, g4->pz,
            gamma1, gamma2, delta);

    // consume arrays
    double sum = 0.0;
    for(unsigned i=0; i<=suml; i++) {
        for(unsigned j=0; j<=summ; j++) {
            for(unsigned k=0; k<=sumn; k++) {
                sum += bx[i]*by[j]*bz[k]*fgamma(i+j+k,0.25*rpq2/delta);
            }
        }
    }

    // delete arrays
    delete[] bx;
    delete[] by;
    delete[] bz;

    return 2.0 * pow(M_PI, 2.5)/(gamma1 * gamma2 * sqrt(gamma1 + gamma2)) *
        exp(-g1->alpha * g2->alpha*rab2/gamma1) *
        exp(-g3->alpha * g4->alpha*rcd2/gamma2) * sum;
}

void build_B_array(double *arr, int l1,int l2,int l3,int l4,
        double p, double a, double b, double q, double c, double d,
        double g1, double g2, double delta) {

    for(int i1=0; i1<l1+l2+1; i1++) {
        for(int i2=0; i2<l3+l4+1; i2++) {
            for(int r1=0; r1 < i1/2+1; r1++) {
                for(int r2=0; r2 < i2/2+1; r2++) {
                    for(int u=0; u<(i1+i2)/2-r1-r2+1; u++) {
                        int i = i1+i2-2*(r1+r2)-u;
                        arr[i] += B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                          p,a,b,q,c,d,g1,g2,delta);
                    }
                }
            }
        }
    }
}

double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2, int l3, int l4,
              double px, double ax, double bx, double qx, double cx, double dx, double gamma1,
              double gamma2, double delta) {
    return fB(i1,l1,l2,px,ax,bx,r1,gamma1)*
        pow(-1,i2) * fB(i2,l3,l4,qx,cx,dx,r2,gamma2)*
        pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u)*
        mypow(qx-px,i1+i2-2*(r1+r2)-2*u)/
        mypow(delta,i1+i2-2*(r1+r2)-u);
}

double fB(int i, int l1, int l2, double p, double a, double b, int r, double g) {
    return binomial_prefactor(i, l1, l2, p-a, p-b) * B0(i, r, g);
}

double B0(int i, int r, double g) {
    return fact_ratio2(i,r) * pow(4*g,r-i);
}

double fact_ratio2(int a, int b) {
    return fact(a) / fact(b) / fact(a - 2*b);
}

unsigned teindex(unsigned i, unsigned j, unsigned k, unsigned l) {
    if(i < j) {
        swap(i,j);
    }
    if(k < l) {
        swap(k,l);
    }

    unsigned ij = i * (i + 1) / 2 + j;
    unsigned kl = k * (k + 1) / 2 + l;

    if(ij < kl) {
        swap(ij,kl);
    }

    return ij * (ij + 1) / 2 + kl;
}
