#include "hf.h"

HF::HF(const Molecule *_mol) : mol(_mol) {
    this->build_matrices();
    this->alpha = 0.7;
    this->nucrep = this->mol->calculate_nuclear_repulsion();
}

void HF::build_matrices() {
    const unsigned nr_atoms = this->mol->get_nr_atoms();
    const unsigned nr_bfs = this->mol->get_nr_bfs();
    int i,j,k; // iterators
    const CGF *cgfs = this->mol->get_cgfs();

    // calculate number of electrons
    this->nrelec = this->mol->get_nr_elec();

    // construct time objects
    clock_t tic, toc;

    // build overlap matrix
    tic = clock();
    this->S = Matrix(nr_bfs, nr_bfs);
    for(i=0; i<nr_bfs; ++i) {
        S(i,i) = cgf_overlap(cgfs[i], cgfs[i]);
        for(j=i+1; j<nr_bfs; ++j) {
            S(i,j) = S(j,i) = cgf_overlap(cgfs[i], cgfs[j]);
        }
    }
    toc = clock();
    printf("Constructed overlap matrix in %4.2f s\n",
        (double)(toc - tic) / CLOCKS_PER_SEC);

    // build kinetic energy matrix
    tic = clock();
    this->T = Matrix(nr_bfs, nr_bfs);
    for(i=0; i<nr_bfs; ++i) {
        for(j=i; j<nr_bfs; ++j) {
            T(i,j) = T(j,i) = cgf_kinetic(cgfs[i], cgfs[j]);
        }
    }
    toc = clock();
    printf("Constructed kinetic energy matrix in %4.2f s\n",
        (double)(toc - tic) / CLOCKS_PER_SEC);

    // build nuclear energy matrix
    tic = clock();
    this->V = Matrix(nr_bfs, nr_bfs); // auto-initializes to zero
    for(k=0; k<nr_atoms; ++k) {
        const double charge = this->mol->get_atom_charge(k);
        const double cx = this->mol->get_atom_coord(k,0);
        const double cy = this->mol->get_atom_coord(k,1);
        const double cz = this->mol->get_atom_coord(k,2);

        for(i=0; i<nr_bfs; ++i) {
            for(j=i; j<nr_bfs; ++j) {
                const double v = charge * cgf_nuclear(cgfs[i], cgfs[j],
                                                      cx, cy, cz);
                if(i != j) { // avoid double-setting
                    this->V(i,j) += v;
                }
                this->V(j,i) += v;
            }
        }
    }
    toc = clock();
    printf("Constructed nuclear attraction matrix in %4.2f s\n",
        (double)(toc - tic) / CLOCKS_PER_SEC);

    // construct two-electron integrals
    tic = clock();
    this->teints = build_two_electron_integrals(cgfs, nr_bfs);
    toc = clock();
    printf("Constructed two electron integrals in %4.2f s\n",
	(double)(toc - tic) / CLOCKS_PER_SEC);

    // perform eigenvalue decomposition
    tic = clock();
    Matrix Vv = this->S.zeros_like();
    Matrix D = this->S.zeros_like();
    eigh(this->S, &Vv, &D);

    // and create transformation matrix
    for(i=0; i<nr_bfs; ++i) {
        D(i,i) = 1.0 / sqrt(D(i,i));
    }
    this->X = Vv * D;
    toc = clock();
    printf("Constructed transformation matrix in %4.2f s\n",
        (double)(toc - tic) / CLOCKS_PER_SEC);

    // calculate core Hamiltonian matrix
    this->H = T + V;

    // build auxiliary matrices
    this->P = Matrix(nr_bfs, nr_bfs);
    this->G = Matrix(nr_bfs, nr_bfs);
}

void HF::run() {
    int i,j,k,l; // regular iterator
    int itr = 0; // SCF iterator
    const unsigned nr_bfs = this->mol->get_nr_bfs();
    double energy;
    double enold = 0.0;
    double ediff = 1e6;  // some large number
    clock_t tic, toc;    // keep track of computaton time

    // start self-consistent field procedure
    while(ediff > 1e-6 && itr < 100) {
        tic = clock();
        ++itr;  // increment SCF counter
        for(i=0; i<nr_bfs; ++i) {
            for(j=0; j<nr_bfs; ++j) {
                this->G(i,j) = 0.0;
                for(k=0; k<nr_bfs; ++k) {
                    for(l=0; l<nr_bfs; ++l) {
                        const unsigned index1 = teindex(i,j,k,l);
                        const unsigned index2 = teindex(i,k,l,j);
                        this->G(i,j) += this->P(k,l) * (this->teints[index1] -
                                               0.5 * this->teints[index2]);
                    }
                }
            }
        }

        // build Fock matrix in regular basis
        this->F = this->H + this->G;

        // perform basis transformation on Fock matrix
        Matrix Fp = this->X.transpose() * this->F * this->X;

        // perform matrix diagonalization
	    Matrix Cp = Fp.zeros_like();
	    Matrix D = Fp.zeros_like();
        eigh(Fp, &Cp, &D);

        // construct coefficient matrix
        this->C = this->X * Cp;

        // store energy matrix
        this->E = D;

        // calculate energy
        Matrix M = this->H + this->F;
        energy = 0.0; // reset energy
        for(i=0; i<nr_bfs; ++i) {
            for(j=0; j<nr_bfs; ++j) {
                energy += P(j,i) * M(i,j);
            }
        }
        energy = 0.5 * energy + this->nucrep; // add nuclear repulsion term

        // construct new density matrix
        Matrix Pnew = Matrix(nr_bfs, nr_bfs); // zero matrix
        for(i=0; i<nr_bfs; ++i) {
            for(j=0; j<nr_bfs; ++j) {
                for(k=0; k<this->nrelec/2; ++k) {
                    Pnew(i,j) += 2.0 * C(i,k) * C(j,k);
                }
            }
        }

        // perform matrix mixing
        for(i=0; i<nr_bfs; ++i) {
            for(j=0; j<nr_bfs; ++j) {
               this->P(i,j) = this->alpha * Pnew(i,j) +
                                (1.0 - alpha) * this->P(i,j);
            }
        }

        // determine energy difference
        ediff = fabs(energy - enold);
        enold = energy;
        toc = clock();
        double t = (double)(toc - tic) / CLOCKS_PER_SEC;
        printf("Iteration %02i | E = %+6.4f Ht | DE = %+6.4f Ht | t = %4.2f s\n",
            itr, energy, ediff, t);
        this->etotal = energy;
    }

    cout << "Calculation complete!" << endl;
}

HF::~HF() {
    free(this->teints);
}