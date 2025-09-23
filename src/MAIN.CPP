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

int main(int argc, char *argv[]) {
    clrscr();

    if(argc < 2) {
        printf("Usage: %s input.mol\n");
        return -1;
    }

    Molecule mol(argv[1]);
    clock_t tic, toc; // keep track of time

    printf("========================\n");
    printf("  Welcome to HF4MSDOS!\n");
    printf("  version: %i.%i.%i\n",
           VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
    printf("========================\n");

    // build and assign basis set
    tic = clock();
    printf("Starting simulation\n");
    mol.print();
    BasisBuilder bb;
    Basis basis = bb.build_basis(mol);
    mol.assign_basis(basis.size(), basis.get_cgfs());
    HF hf(&mol);
    hf.run();
    toc = clock();
    double calctime = (double)(toc - tic) / CLOCKS_PER_SEC;

    char c;
    printf("Simulation complete (%.2f seconds)\n", calctime);
    while(1) {
        printf("Which output would you like to view?\n");
        printf("0 - Print molecular coordinates\n");
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
        while(!((c>='0' && c<='8') || c == 'q')) {
            c = getch();
        }
        switch(c) {
            case '0':
                mol.print();
                wait_for_key();
            break;
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