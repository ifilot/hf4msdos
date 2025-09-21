#include "mymath.h"

double fact(unsigned x) {
    static const double lookup[11] = {
        1.0,
        1.0,
        2.0,
        6.0,
        24.0,
        120.0,
        720.0,
        5040.0,
        40320.0,
        362880.0,
        3628800.0
    };

    if(x <= 10) {
        return lookup[x];
    } else {
        double acc = lookup[10];
        for(unsigned k = 11u; k <= x; ++k) {
            acc *= (double)k;
        }
        return acc;
    }
}

double fact2(int x) {
    static const double lookup[11] = {
        1.0,
        1.0,
        2.0,
        3.0,
        8.0,
        15.0,
        48.0,
        105.0,
        384.0,
        945.0,
        3840.0
    };

    if(x < -1) {
        cerr << "Invalid double factorial!" << endl;
        exit(-1);
    }

    if(x == -1) {
        return 1.0;
    }

    if(x <= 11) {
        return lookup[x];
    } else {
        double acc = ((x & 1u) == 0u) ? lookup[10] : lookup[9];
        unsigned k0 = ((x & 1u) == 0u) ? 12u : 11u;
        for(unsigned k = k0; k<=x; k += 2u) {
            acc *= (double)k;
        }
        return acc;
    }
}

double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb) {
    double sum = 0;

    for(int t=0; t<(s+1); ++t) {
        if((s-ia <= t) && (t <= ib)) {
            // the code below prevents Borland C++ from
            // barking POW DOMAIN ERROR
            double t0 = 0.0;
            if(ia-s+t == 0) {
                t0 = 1.0;   // 0^0 is thus defined as 1.0 here
            } else if(xpa == 0) {
                continue;   // result is given, so we ignore this term
            } else {
                t0 = pow(xpa, ia-s+t);
            }

            double t1 = 0.0;
            if(ib-t == 0) {
                t1 = 1.0;
            } else if(xpb == 0) {
                continue;
            } else {
                t1 = pow(xpb, ib-t);
            }


            sum += binomial(ia, s-t) * binomial(ib, t) * t0 * t1;
        }
    }

    return sum;
}

int binomial(int a, int b) {
    return fact(a) / (fact(b) * fact(a-b));
}

double mypow(double base, double exp) {
    if(exp == 0.0) {
        return 1.0;
    } else if(base == 0.0) {
        return 0.0;
    } else {
        return pow(base, exp);
    }
}

void swap(unsigned &a, unsigned &b) {
    unsigned c = a;
    a = b;
    b = c;
}