#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;


// global variables

#define fichierLecture "parameters.csv"     // names of input

// definition of structure "chr" representing a chromosome:
// "S" is the S-allele
// "sel" is a vector containing the positions of deleterious alleles along the chromosome
// "let" is a vector containing the positions of lethal alleles along the chromosome

// a chromosome
struct chr
{
   int S; //S-locus
   vector<double> sel; // selected loci: vector of deleterious mutations
   vector<double> let; // vector of lethal mutations
   vector<double> neutral_S; // vector of neutral allele linked to S locus
   vector<double> neutral_ind; // vector of neutral allele independent of S locus

};

// results of the recursion function for the main function
struct result
{
       double frSC; // frequency of SC
       double delta; // inbreeding depression
       };

// for some mesure of genetic structure
struct resultats_FST
{
       double FST;
};

// Prototypes of functions

void ouvrirFichierE(); // open input file
void fermerFichierE(); // close input file

// read input file
bool lireFichier(int &Numr, int &Nr,int &pr, int &max_alleles, double &dpr, double &dgr, double &ar, double &sr, double &hr, double &slr, double &hlr, double &thr, double &gammar,
                     double &Uir, int &itr, double &pasir,double &seuilr, double &Usir, double &Uscr, double &Uner, int &mUscr, int &menager, double &Lr,
                     int &NbPrelimAr, int &NbPrelimBr, int &NbGenr, int &pasprr, int &passcr, int &sampleSr, int &nbr_locus_neutrer, int &position_stream, int &etat);

// function of systeme evolution
result recursion(int numsimulv, int iv, int Nv, double av, double sv, double hv, double slv, double hlv, double thv, double gammav, double Uv,
                     double Usiv, double Uscv, double Unev, int mUscv, double Lv, double seuilv, int NbPrelimAv, int NbPrelimBv,
					 int NbGenv, int pasprv, int pasv, int nbr_locus_neutrev, int sampleSv,  double dpv , double dgv , int pv, int menagev, int nbr_allelesv);

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);

double fitness(chr &c1, chr &c2, double wHe, double wHo, double wHe_l, double wHo_l); //fitness 
void rec(chr &res, chr &c1, chr &c2, double sz); // recombination
int menage_del(chr * pop,int taille_pop_chr); // delete fixed deleterious mutations
int menage_let(chr * pop,int taille_pop_chr); // delete fixed lethals mutations
void cntl_c_handler(int bidon); // manualy stop

#endif
