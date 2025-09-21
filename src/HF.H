#ifndef _HF_H
#define _HF_H

#include <time.h>

#include "molecule.h"
#include "matrix.h"
#include "integrals.h"
#include "teints.h"
#include "eig.h"

class HF {

private:
    const Molecule *mol;

    Matrix S; // overlap matrix
    Matrix T; // kinetic energy matrix
    Matrix V; // nuclear attraction matrix
    Matrix H; // core Hamiltonian matrix
    Matrix X; // transformation matrix
    Matrix P; // density matrix
    Matrix G; // auxiliary matrix
    Matrix F; // Fock matrix
    Matrix C; // coefficient matrix
    Matrix E; // energy matrix

    unsigned nrelec; // total number of electrons in molecule
    double* teints;  // array of two-electron integrals
    double alpha;    // mixing constant
    double nucrep;   // nuclear repulsion energy
    double etotal;   // total electronic energy

public:
    HF(const Molecule *_mol);

    inline const Matrix& get_overlap_matrix() const {
        return this->S;
    }

    inline const Matrix& get_kinetic_matrix() const {
        return this->T;
    }

    inline const Matrix& get_nuclear_matrix() const {
        return this->V;
    }

    inline const Matrix& get_fock_matrix() const {
        return this->F;
    }

    inline const Matrix& get_coefficient_matrix() const {
        return this->C;
    }

    inline const Matrix& get_energy_matrix() const {
        return this->E;
    }

    inline double get_total_electronic_energy() const {
        return this->etotal;
    }

    void run();

    ~HF();

private:
    void build_matrices();

};

#endif