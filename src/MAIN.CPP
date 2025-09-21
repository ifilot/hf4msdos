#include <iostream.h>
#include <conio.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "molecule.h"
#include "basisset.h"
#include "cgf.h"
#include "integrals.h"
#include "matrix.h"
#include "teints.h"
#include "eig.h"
#include "hf.h"
#include "version.h"

/*
 * Helper function that asks user to press any key to continue
 */
static void wait_for_key() {
    printf("-- Press any key to continue --");
    getch();
    printf("\r");
}

int main() {
    clrscr();

    /* NOTE: All atom coordinates are in Angstrom but will be
             automatically converted to Bohr units (a.u.) in the
             Molecule class */

    // H2 molecule
    Molecule h2(2);
    h2.set_atom(0, "H", 0.0, 0.0, -0.3679020);
    h2.set_atom(1, "H", 0.0, 0.0,  0.3679020);
    //h2.print();

    // CO molecule
    Molecule co(2);
    co.set_atom(0, "C", 0, 0,  0.5120840);
    co.set_atom(1, "O", 0, 0, -0.6827790);

    // LiH molecule
    Molecule lih(2);
    lih.set_atom(0, "Li", 0.0, 0.0,  0.724572);
    lih.set_atom(1, "H" , 0.0, 0.0, -2.173714);

    Molecule h2o(3);
    h2o.set_atom(0, "O", 0, 0, 0.1368220);
    h2o.set_atom(1, "H", 0, 0.7697660, -0.5472890);
    h2o.set_atom(2, "H", 0, -0.7697660, -0.5472890);

    clock_t tic, toc; // keep track of time

    printf("========================\n");
    printf("  Welcome to HF4MSDOS!\n");
    printf("  version: %i.%i.%i\n",
           VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
    printf("========================\n");
    printf("I come pre-compiled with a number of molecules.\n");
    printf("Which molecule would you like me to simulate?\n");
    printf("1 - H2\n");
    printf("2 - CO\n");
    printf("3 - LiH\n");
    printf("4 - H2O\n");
    printf("Q - QUIT\n");
    printf("Press 1-4 to start simulation or q to quit.\n");

    char c = getch();
    while(!( (c>='1' && c<='4') || c == 'q')) {
        c = getch();
    }

    Molecule *mol;
    switch(c) {
        case '1':
            mol = &h2;
        break;
        case '2':
            mol = &co;
        break;
        case '3':
            mol = &lih;
        break;
        case '4':
            mol = &h2o;
        break;
        case 'q':
            printf("Goodbye!\n");
            return 0;
    }

    // build and assign basis set
    tic = clock();
    printf("\nStarting simulation\n");
    BasisBuilder bb;
    Basis basis = bb.build_basis(*mol);
    mol->assign_basis(basis.size(), basis.get_cgfs());
    HF hf(mol);
    hf.run();
    toc = clock();
    double calctime = (double)(toc - tic) / CLOCKS_PER_SEC;

    printf("Simulation complete (%.2f seconds)\n", calctime);
    while(1) {
        printf("Which output would you like to view?\n");
        printf("1 - Overlap matrix\n");
        printf("2 - Kinetic matrix\n");
        printf("3 - Nuclear matrix\n");
        printf("4 - Fock matrix\n");
        printf("5 - Coefficient matrix\n");
        printf("6 - Energy matrix\n");
        printf("7 - Total electronic energy\n");
        printf("8 - Total calculation time\n");
        printf("q - Quit\n");
        c = getch();
        while(!((c>='1' && c<='8') || c == 'q')) {
            c = getch();
        }
        switch(c) {
            case '1':
                printf("Overlap matrix:\n");
                hf.get_overlap_matrix().print("% 6.4f ");
                wait_for_key();
            break;
            case '2':
                printf("Kinetic matrix:\n");
                hf.get_kinetic_matrix().print("% 6.4f ");
                wait_for_key();
            break;
            case '3':
                printf("Nuclear matrix:\n");
                hf.get_nuclear_matrix().print("% 6.4f ");
                wait_for_key();
            break;
            case '4':
                printf("Fock matrix:\n");
                hf.get_fock_matrix().print("% 6.4f ");
                wait_for_key();
            break;
            case '5':
                printf("Coefficient matrix:\n");
                hf.get_coefficient_matrix().print("% 6.4f ");
                wait_for_key();
            break;
            case '6':
                printf("Energy matrix:\n");
                hf.get_energy_matrix().print("% 6.4f ");
                wait_for_key();
            break;
            case '7':
                printf("Total electronic energy:\n");
                printf("Etotal = %12.8f Ht\n", hf.get_total_electronic_energy());
            break;
            case '8':
                printf("Total calculation time: %.2f seconds\n", calctime);
            break;
            case 'q':
                printf("Goodbye!");
                return 0;
        }
    }
}