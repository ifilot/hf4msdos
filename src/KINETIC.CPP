#include "kinetic.h"

/*
 * Calculate the kinetic integral between two CGF objects
 */
double cgf_kinetic(const CGF& cgf1, const CGF& cgf2) {
    double sum = 0;
    const unsigned char N1 = cgf1.get_number_gtos();
    const unsigned char N2 = cgf2.get_number_gtos();

    for(unsigned char i=0; i<N1; ++i) {
	const struct GTO *g1 = cgf1.get_gto_ptr(i);
	const double c1 = cgf1.get_coeff(i);
	for(unsigned char j=0; j<N2; ++j) {
	    const struct GTO *g2 = cgf2.get_gto_ptr(j);
	    const double c2 = cgf2.get_coeff(j);
	    double term = c1 * c2 * gto_kinetic(g1, g2);
            sum += term;
        }
    }

    return sum;
}

/*
 * Calculate kinetic energy integral of two GTOs
 */
double gto_kinetic(const struct GTO *g1, const struct GTO *g2) {
    const double t0 = g2->alpha * (2 * (g2->l + g2->m + g2->n) + 3) *
                      gto_overlap(g1, g2);

    const double t1 = -2 * g2->alpha * g2->alpha * (
        gto_overlap_modcoeff(g1, g2, g1->l, g1->m, g1->n, g2->l+2, g2->m, g2->n) +
        gto_overlap_modcoeff(g1, g2, g1->l, g1->m, g1->n, g2->l, g2->m+2, g2->n) +
        gto_overlap_modcoeff(g1, g2, g1->l, g1->m, g1->n, g2->l, g2->m, g2->n+2)
    );

    const double t2 = -0.5 * (g2->l * (g2->l - 1) *
        gto_overlap_modcoeff(g1, g2, g1->l, g1->m, g1->n, g2->l-2, g2->m * (g2->m - 1), g2->n) +
        gto_overlap_modcoeff(g1, g2, g1->l, g1->m, g1->n, g2->l, g2->m-2, g2->n*(g2->n-1)) +
        gto_overlap_modcoeff(g1, g2, g1->l, g1->m, g1->n, g2->l, g2->m, g2->n-2)
    );

    return t0 + t1 + t2;
}