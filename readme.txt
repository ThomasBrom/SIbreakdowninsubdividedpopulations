The general purpose of the program is to investigate the maintenance of a gametophytic self-incompatibility system in front of self-compatible mutants.
The program do evolutionary simulations for different value of deleterious (or lethal) mutations that lead to inbreeding depression.

recursion.cpp is the function for each simulation, that receive a set of parameters to test the invasion (or not) by the SC mutants

main.cpp use fichier.cpp to read the parameters.csv file and to launche "recursion".
main can investigate different value of U sequencialy or search for a critical value of U: the lowest value for which SI system is maintained

Mercennes Twister.h and ranbin.cpp are for random numbers

fichiers.cpp is to open, close, read parameters file.

SelRec.cpp contain mainly functions for fitness calculation and for recombination

depression.h contain function and structurs declaration

contact: Thomas Brom (thomas.brom@univ-lille.fr or researchgate.net/profile/Thomas_Brom or find me on web)
