#include "../header/Interface.h"

/**
 * @author Sergey Dyachkov.
 * @date 23.12.2013.
 * @mainpage Thomas Fermi Model with quantum and exchange corrections.
 *
 * @section short_description Short description.
 * 
 * This program is being developed for calculation of thermodynamic functions of electrons
 * followed from Thomas-Fermi model with quantum and exchange corrections. The Interface 
 * was developed to read file as command line argument and calculate the values of functions
 * that are indicated in the file.
 */

int main(int argc, char** argv) {

    Interface i;
    i.Read(argv[1]);
    i.Print();
    return 0;
}

