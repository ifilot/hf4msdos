#ifndef _CGF_H
#define _CGF_H

#include <alloc.h>
#include <iostream.h>
#include <stdlib.h>
#include <conio.h>
#include <stdio.h>

#include "gto.h"

/*
 * Collects multiple GTOs to form a Contracted Gaussian Functional
 *
 * This class is relatively flexible, it allows for GTOs located
 * at different centers (relevant for SALCs). GTO coefficients
 * include the normalization constant of the GTO, i.e., they are
 * the product of the CGF expansion coefficient **and** the GTO
 * normalization constant.
 */
class CGF {
private:
    unsigned char N;         /* number of GTO in CGF */
    GTO *gto_list;           /* pointer to first element of GTO */
    double *gto_coeff;       /* pointer to first element of GTO coefficients */

public:
    CGF();
    /*
     * Create a new CGF
     *
     * @param N : number of GTOs (fixed)
     */
    CGF(unsigned _N);

    void set_ngtos(unsigned _N);

    /*
     * Set a GTO in the CGF
     *
     * @param idx : GTO index
     * @param coeff : GTO coefficient in CGF
     * @param n : power in x
     * @param l : power in y
     * @param m : power in z
     * @param alpha : exponent
     * @param px : x-position of GTO
     * @param py : y-position of GTO
     * @param pz : z-position of GTO
     */
    void set_gto(unsigned idx,
                 double coeff,
                 int n,
                 int l,
                 int m,
                 double alpha,
                 double px,
                 double py,
                 double pz);

    /*
     * Get the number of GTOs stored in CGF
     */
    inline const unsigned char get_number_gtos() const {
        return this->N;
    }

    /*
     * Get a constant pointer to a GTO
     */
    inline const struct GTO* get_gto_ptr(unsigned idx) const {
        return &this->gto_list[idx];
    }

    /*
     * Get the expansion coefficient for the GTO, including the
     * GTOs normalization constant.
     */
    inline const double get_coeff(unsigned idx) const {
        return this->gto_coeff[idx];
    }

    void print() const;
};

#endif