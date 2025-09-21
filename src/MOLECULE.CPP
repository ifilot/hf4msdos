#include "molecule.h"

// this array is used to convert element names to element
// numbers (used for indexing and for determining charges)
static const char *elnames[] = {
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne"
};

Molecule::Molecule(unsigned _nr_atoms) : nr_atoms(_nr_atoms) {
    if(this->nr_atoms > 0) {
        this->atom_pos = (double*)calloc(this->nr_atoms*3, sizeof(double));
        this->atom_el = (unsigned char*)calloc(this->nr_atoms, sizeof(unsigned char));
    }
}

void Molecule::set_atom(unsigned elid, const char *el, double px, double py, double pz) {
    if(elid >= this->nr_atoms) {
        cerr << "Invalid element id encountered!" << endl;
        return;
    }
    this->atom_pos[elid * 3] = px * ANG2BOHR;
    this->atom_pos[elid * 3 + 1] = py * ANG2BOHR;
    this->atom_pos[elid * 3 + 2] = pz * ANG2BOHR;

    for(int i=0; i<10; i++) {
        if(strcmp(el,elnames[i]) == 0) {
            this->atom_el[elid] = i+1;
            return;
        }
    }

    cerr << "Could not identify atom!" << endl;
}

unsigned Molecule::get_nr_elec() const {
    unsigned nrelec = 0;
    for(int i=0; i<this->nr_atoms; ++i) {
        nrelec += this->get_atom_charge(i);
    }

    return nrelec;
}

Molecule::~Molecule() {
    if(this->nr_atoms > 0) {
        free(this->atom_pos);
        free(this->atom_el);
    }
}

void Molecule::print() {
    for(int i=0; i<this->nr_atoms; ++i) {
        printf("%2s  %6.4f  %6.4f  %6.4f\n",
               elnames[this->atom_el[i]-1],
               this->atom_pos[i*3+0],
               this->atom_pos[i*3+1],
               this->atom_pos[i*3+2]);
    }
}

double Molecule::get_atom_distance(unsigned i, unsigned j) const {
    const double ax = this->get_atom_coord(i,0);
    const double ay = this->get_atom_coord(i,1);
    const double az = this->get_atom_coord(i,2);

    const double bx = this->get_atom_coord(j,0);
    const double by = this->get_atom_coord(j,1);
    const double bz = this->get_atom_coord(j,2);

    const double dx = (ax - bx);
    const double dy = (ay - by);
    const double dz = (az - bz);

    return sqrt(dx*dx + dy*dy + dz*dz);
}

double Molecule::calculate_nuclear_repulsion() const {
    double nucrep = 0.0;

    for(int i=0; i<this->nr_atoms; ++i) {
        const double c1 = this->get_atom_charge(i);
        for(int j=i+1; j<this->nr_atoms; ++j) {
            const double c2 = this->get_atom_charge(j);
            const double d = get_atom_distance(i,j);
            nucrep += 1.0 / d * c1 * c2;
        }
    }

    return nucrep;
}

void Molecule::assign_basis(unsigned sz, const CGF* _cgfs) {
    this->nr_bfs = sz;
    this->cgfs = _cgfs;
}