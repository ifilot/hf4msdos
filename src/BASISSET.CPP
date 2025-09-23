#include "basisset.h"

/*
 * Specifies which shells an atom has
 */
static const unsigned char shells[] = {
    SHELL_1S,                       // H
    SHELL_1S,                       // He
    SHELL_1S | SHELL_2S | SHELL_2P, // Li
    SHELL_1S | SHELL_2S | SHELL_2P, // Be
    SHELL_1S | SHELL_2S | SHELL_2P, // B
    SHELL_1S | SHELL_2S | SHELL_2P, // C
    SHELL_1S | SHELL_2S | SHELL_2P, // N
    SHELL_1S | SHELL_2S | SHELL_2P, // O
    SHELL_1S | SHELL_2S | SHELL_2P, // F
    SHELL_1S | SHELL_2S | SHELL_2P  // Ne
};

/*
 * Basis set coefficients for STO-3g stored as 6-tuples. Should
 * be read as 3 x (alpha,coeff) pairs.
 */
static double param_1s[][6] = {
    {  3.425251, 0.154329,  0.623914, 0.535328, 0.168855, 0.444635}, // H
    {  6.362421, 0.154329,  1.158923, 0.535328, 0.313650, 0.444635}, // He
    { 16.119575, 0.154329,  2.936201, 0.535328, 0.794651, 0.444635}, // Li
    {000.000000, 0.000000, 00.000000, 0.000000, 0.000000, 0.000000}, // Be
    {000.000000, 0.000000, 00.000000, 0.000000, 0.000000, 0.000000}, // B
    { 71.616837, 0.154329, 13.045096, 0.535328, 3.530512, 0.444635}, // C
    { 99.106169, 0.154329, 18.052312, 0.535328, 4.885660, 0.444635}, // N
    {130.70932 , 0.154329, 23.808861, 0.535328, 6.443608, 0.444635}, // O
    {000.000000, 0.000000, 00.000000, 0.000000, 0.000000, 0.000000}, // F
    {000.000000, 0.000000, 00.000000, 0.000000, 0.000000, 0.000000}  // Ne
};

static double param_2s[][6] = {
    {3.425251, 0.154329, 0.623914, 0.535328, 0.168855, 0.444635}, // H
    {6.362421, 0.154329, 1.158923, 0.535328, 0.313650, 0.444635}, // He
    {0.63629, -0.099967, 0.147860, 0.399513, 0.048089, 0.700115}, // Li
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // Be
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // B
    {2.941249,-0.099967, 0.683483, 0.399513, 0.222290, 0.700115}, // C
    {3.780456,-0.099967, 0.878497, 0.399513, 0.285714, 0.700115}, // N
    {5.033151,-0.099967, 1.169596, 0.399513, 0.380389, 0.700115}, // O
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // F
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}  // Ne
};

static double param_2p[][6] = {
    {3.425251, 0.154329, 0.623914, 0.535328, 0.168855, 0.444635}, // H
    {6.362421, 0.154329, 1.158923, 0.535328, 0.313650, 0.444635}, // He
    {0.636290, 0.155916, 0.147860, 0.607684, 0.048089, 0.391957}, // Li
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // Be
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // B
    {2.941249, 0.155916, 0.683483, 0.607684, 0.222290, 0.391957}, // C
    {3.780456, 0.115916, 0.878497, 0.607684, 0.285714, 0.391957}, // N
    {5.033151, 0.155916, 1.169596, 0.607684, 0.380389, 0.391957}, // O
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, // F
    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}  // Ne
};

BasisBuilder::BasisBuilder() {}

Basis BasisBuilder::build_basis(const Molecule& mol) {
    const unsigned natoms = mol.get_nr_atoms();
    unsigned norbs = 0;
    int i,j,k;

    // calculate total number of CGFs
    for(i=0; i<natoms; ++i) {
        const unsigned elid = mol.get_atom_charge(i);
        if(shells[elid] & SHELL_1S) {
            ++norbs;
        }
        if(shells[elid] & SHELL_2S) {
            ++norbs;
        }
        if(shells[elid] & SHELL_2P) {
            norbs += 3;
        }
    }

    // build CGFs
    CGF *cgfs = new CGF[norbs];
    unsigned orbctr = 0;
    for(i=0; i<natoms; ++i) { // loop over atoms
        const unsigned elid = mol.get_atom_charge(i)-1; // note: atom id is charge - 1

        // 1S CGFs
        if(shells[elid] & SHELL_1S) {
            cgfs[orbctr] = CGF(3);
            for(j=0; j<3; ++j) {  // loop over GTOs
                 cgfs[orbctr].set_gto(j, param_1s[elid][j*2+1], // coeff
                     0, 0, 0, param_1s[elid][j*2+0], // alpha
                     mol.get_atom_coord(i,0),
                     mol.get_atom_coord(i,1),
                     mol.get_atom_coord(i,2));
            }
            ++orbctr;
        }
        // 2S CGFs
        if(shells[elid] & SHELL_2S) {
            cgfs[orbctr] = CGF(3);
            for(j=0; j<3; ++j) {  // loop over GTOs
                 cgfs[orbctr].set_gto(j, param_2s[elid][j*2+1], // coeff
                     0, 0, 0, param_2s[elid][j*2+0], // alpha
                     mol.get_atom_coord(i,0),
                     mol.get_atom_coord(i,1),
                     mol.get_atom_coord(i,2));
            }
            ++orbctr;
        }

        // 2P CGFs
        if(shells[elid] & SHELL_2P) {
            for(k=0; k<3; ++k) { // loop over l values
                int nn[3] = {0,0,0};
                nn[k] += 1;
                cgfs[orbctr] = CGF(3);
                for(j=0; j<3; ++j) { // loop over GTOs
                     cgfs[orbctr].set_gto(j, param_2p[elid][j*2+1], // coeff
                         nn[0], nn[1], nn[2], param_2p[elid][j*2+0],// alpha
                         mol.get_atom_coord(i,0),
                         mol.get_atom_coord(i,1),
                         mol.get_atom_coord(i,2));
                }
                ++orbctr;
            }
        }
    }
    return Basis(orbctr, cgfs);
}

Basis::Basis(const Basis& _other) {
    this->N = _other.N;
    this->cgfs = _other.cgfs;
}

Basis::Basis(unsigned _N, CGF* _cgfs) : N(_N), cgfs(_cgfs) {}

Basis::~Basis() {
    delete[] cgfs;
}