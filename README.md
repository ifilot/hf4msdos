# Hartree-Fock for MS-DOS

This repository contains a localized-orbital spin-restricted Hartree-Fock
calculator. The code has been built using the Borland C++ compiler using a
compact memory model and has been tested to work on an IBM 5150 with 256 KiB of
memory.

## Compilation instructions

Install [Borland's Turbo C++](https://winworldpc.com/product/borland-c/30) on
a compatible MS-DOS system. It is recommended to use a 80386 or better with
2 MiB of memory. Open `HF4MSDOS.PRJ` and use the `Build All` function to
compile the program.

## Compilation settings

Depending on your desired deployment platform, you might want to adjust the
compiler settings to produce code for 8086, 80186, or 80286 processors and/or
whether to use a software math library or to compile for a 8087, 80187 or 80287
floating point co-processor unit.