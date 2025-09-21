#include "overlap.h"

/*
 * Calculate the overlap integral between two CGF objects
 */
double cgf_overlap(const CGF& cgf1, const CGF& cgf2) {
    double sum = 0;
    const unsigned char N1 = cgf1.get_number_gtos();
    const unsigned char N2 = cgf2.get_number_gtos();

    for(unsigned char i=0; i<N1; ++i) {
        const struct GTO *g1 = cgf1.get_gto_ptr(i);
        const double c1 = cgf1.get_coeff(i);
        for(unsigned char j=0; j<N2; ++j) {
            const struct GTO *g2 = cgf2.get_gto_ptr(j);
            const double c2 = cgf2.get_coeff(j);
            double term = c1 * c2 * gto_overlap(g1, g2);
            sum += term;
	    }
    }

    return sum;
}

/*
 * Calculate overlap integral of two GTOs
 */
double gto_overlap(const struct GTO *g1, const struct GTO *g2) {
    const double dx = g1->px - g2->px;
    const double dy = g1->py - g2->py;
    const double dz = g1->pz - g2->pz;

    const double rab2 = (dx*dx) + (dy*dy) + (dz*dz);
    const double gamma = g1->alpha + g2->alpha;

    const double gx = (g1->alpha * g1->px + g2->alpha * g2->px) / gamma;
    const double gy = (g1->alpha * g1->py + g2->alpha * g2->py) / gamma;
    const double gz = (g1->alpha * g1->pz + g2->alpha * g2->pz) / gamma;

    const double pre = pow(M_PI / gamma, 1.5) * exp(-g1->alpha * g2->alpha * rab2 / gamma);
    const double wx = gto_overlap_1d(g1->l, g2->l, gx - g1->px, gx - g2->px, gamma);
    const double wy = gto_overlap_1d(g1->m, g2->m, gy - g1->py, gy - g2->py, gamma);
    const double wz = gto_overlap_1d(g1->n, g2->n, gz - g1->pz, gz - g2->pz, gamma);

    return pre * wx * wy * wz;
}


/*
 * Calculate overlap integral of two GTOs with modified Cartesian powers.
 * This function is used as a subroutine in the kinetic overlap integral
 */
double gto_overlap_modcoeff(const struct GTO *g1, const struct GTO *g2,
                            int l1, int m1, int n1,
                            int l2, int m2, int n2) {
    const double dx = g1->px - g2->px;
    const double dy = g1->py - g2->py;
    const double dz = g1->pz - g2->pz;

    const double rab2 = (dx*dx) + (dy*dy) + (dz*dz);
    const double gamma = g1->alpha + g2->alpha;

    const double gx = (g1->alpha * g1->px + g2->alpha * g2->px) / gamma;
    const double gy = (g1->alpha * g1->py + g2->alpha * g2->py) / gamma;
    const double gz = (g1->alpha * g1->pz + g2->alpha * g2->pz) / gamma;

    const double pre = pow(M_PI / gamma, 1.5) * exp(-g1->alpha * g2->alpha * rab2 / gamma);
    const double wx = gto_overlap_1d(l1, l2, gx - g1->px, gx - g2->px, gamma);
    const double wy = gto_overlap_1d(m1, m2, gy - g1->py, gy - g2->py, gamma);
    const double wz = gto_overlap_1d(n1, n2, gz - g1->pz, gz - g2->pz, gamma);

    return pre * wx * wy * wz;
}

/*
 * Calculate overlap integral for one Cartesian component; used in
 * GTO overlap
 */
double gto_overlap_1d(int l1, int l2, double x1, double x2, double gamma) {
    double sum = 0;

    for(int i=0; i < (1 + floor(0.5 * (l1 + l2))); ++i) {
        sum += binomial_prefactor(2*i, l1, l2, x1, x2) *
               fact2(2*i-1) / pow(2 * gamma, i);
    }

    return sum;
}