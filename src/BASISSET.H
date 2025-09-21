#ifndef _BASISSET_H
#define _BASISSET_H

#include "molecule.h"

#define SHELL_1S (1 << 0)
#define SHELL_2S (1 << 1)
#define SHELL_2P (1 << 2)

class Basis {
private:
    unsigned N; // size of the basis
    CGF* cgfs;  // pointer to CGF array
public:
    Basis(const Basis& _other);
    Basis(unsigned _N, CGF* _cgfs);
    ~Basis();

    inline unsigned size() const { return this->N; }
    inline const CGF* get_cgfs() const { return this->cgfs; }

};

class BasisBuilder {
private:

public:
    BasisBuilder();

    struct Basis build_basis(const Molecule& mol);

private:
};

#endif