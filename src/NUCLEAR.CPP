#include "nuclear.h"

/*
 * Calculate the nuclear integral between two CGF objects
 */
double cgf_nuclear(const CGF& cgf1, const CGF& cgf2,
                   double cx, double cy, double cz) {
    double sum = 0;
    const unsigned char N1 = cgf1.get_number_gtos();
    const unsigned char N2 = cgf2.get_number_gtos();

    for(unsigned char i=0; i<N1; ++i) {
	const struct GTO *g1 = cgf1.get_gto_ptr(i);
	const double c1 = cgf1.get_coeff(i);
	for(unsigned char j=0; j<N2; ++j) {
	    const struct GTO *g2 = cgf2.get_gto_ptr(j);
	    const double c2 = cgf2.get_coeff(j);
	    double term = c1 * c2 * gto_nuclear(g1, g2, cx, cy, cz);
            sum += term;
        }
    }

    return sum;
}

/*
 * Calculate nuclear integral for two GTOs
 */
double gto_nuclear(const struct GTO *g1, const struct GTO *g2,
                   double nx, double ny, double nz) {
    const double gamma = g1->alpha + g2->alpha;

    const double dx = g1->px - g2->px;
    const double dy = g1->py - g2->py;
    const double dz = g1->pz - g2->pz;
    const double rab2 = (dx * dx + dy * dy + dz * dz);

    // calculate Gaussian product center
    const double gx = (g1->alpha * g1->px + g2->alpha * g2->px) / gamma;
    const double gy = (g1->alpha * g1->py + g2->alpha * g2->py) / gamma;
    const double gz = (g1->alpha * g1->pz + g2->alpha * g2->pz) / gamma;

    const double cx = nx - gx;
    const double cy = ny - gy;
    const double cz = nz - gz;
    const double rcp2 = (cx * cx + cy * cy + cz * cz);

    // create arrays and initialize them to zero
    unsigned il = g1->l + g2->l + 1;
    unsigned im = g1->m + g2->m + 1;
    unsigned in = g1->n + g2->n + 1;
    double *ax = new double[il];
    memset(ax, 0.0, sizeof(double) * il);
    double *ay = new double[im];
    memset(ay, 0.0, sizeof(double) * im);
    double *az = new double[in];
    memset(az, 0.0, sizeof(double) * in);

    // populate arrays
    build_A_array(ax, il, g1->l, g2->l, gx-g1->px, gx-g2->px, gx-nx, gamma);
    build_A_array(ay, im, g1->m, g2->m, gy-g1->py, gy-g2->py, gy-ny, gamma);
    build_A_array(az, in, g1->n, g2->n, gz-g1->pz, gz-g2->pz, gz-nz, gamma);

    // consume arrays
    double sum = 0.0;
    for(int i=0; i<=g1->l+g2->l;++i) {
        for(int j=0; j<=g1->m+g2->m;++j) {
            for(int k=0; k<=g1->n+g2->n; ++k) {
                sum += ax[i] * ay[j] * az[k] * fgamma(i+j+k, rcp2*gamma);
            }
        }

    }

    // free the dynamically allocated arrays
    delete[] ax;
    delete[] ay;
    delete[] az;

    return -2.0 * M_PI / gamma * exp(-g1->alpha * g2->alpha * rab2 / gamma) * sum;
}

/*
 * Build A array for nuclear integral. Allocates memory to populate result
 * vector. Note that allocated memory needs to be free'd later.
 */
void build_A_array(double *arr, unsigned imax,
                   int l1, int l2,
                   double pa, double pb, double cp, double g) {
    for(int i=0; i<imax; ++i) {
        for(int r=0; r<=i/2; ++r) {
            for(int u=0; u<=(i-2*r)/2; ++u) {
                int iI = i - 2 * r - u;
                arr[iI] += a_term(i, r, u, l1, l2, pa, pb, cp, g);
            }
        }
    }
}

double a_term(int i, int r, int u, int l1, int l2,
              double pax, double pbx, double cpx, double gamma) {

    return pow(-1,i) * binomial_prefactor(i,l1,l2,pax,pbx) *
           pow(-1,u) * fact(i)*mypow(cpx,i-2*r-2*u) *
           pow(0.25/gamma,r+u)/fact(r)/fact(u)/fact(i-2*r-2*u);
}