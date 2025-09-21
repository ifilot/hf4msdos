#ifndef _MOLECULE_H
#define _MOLECULE_H

#include <stdlib.h>
#include <iostream.h>
#include <string.h>
#include <stdio.h>

#include "cgf.h"

#define ANG2BOHR 1.88973

class Molecule {

private:
    unsigned nr_atoms;
    double *atom_pos;
    unsigned char *atom_el;
    unsigned nr_bfs;
    const CGF *cgfs;

public:
    Molecule(unsigned _nr_atoms);

    void set_atom(unsigned elid, const char *el, double px, double py, double pz);

    void print();

    void assign_basis(unsigned sz, const CGF* cgfs);

    inline unsigned get_nr_bfs() const {
        return this->nr_bfs;
    }

    inline unsigned get_nr_atoms() const {
        return this->nr_atoms;
    }

    inline const CGF* get_cgfs() const {
        return this->cgfs;
    }

    inline double get_atom_charge(unsigned elid) const {
        return this->atom_el[elid];
    }

    inline double get_atom_coord(unsigned elid, unsigned i) const {
        return this->atom_pos[elid * 3 + i];
    }

    double get_atom_distance(unsigned i, unsigned j) const;

    double calculate_nuclear_repulsion() const;

    unsigned get_nr_elec() const;

    ~Molecule();

private:

};

#endif