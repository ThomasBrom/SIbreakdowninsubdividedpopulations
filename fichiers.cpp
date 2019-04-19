// Functions to open input and output files,
// read parameter values and write them in output file

#include "depression.h"
#include <iostream>
#include <fstream>
#include <sstream>
#ifdef __unix__
#include <unistd.h>
#else
#include <Windows.h>
#endif
using namespace std;

extern FILE * fichierE;
FILE * fichierblocage;

//Opens input file:
void ouvrirFichierE()
{
 int compteur_tentative=0;
	// 10 attempt
    while (compteur_tentative<10) {

		// look for a "blocking" file
        fichierblocage = fopen("blocage.txt","r");

        if (fichierblocage==NULL) { // if there is no "blocking" file
            ofstream foutblocage("blocage.txt"); // add a blocking file, this process will work on the file !!
            foutblocage.close();

			// open parameters file
            fichierE = fopen(fichierLecture,"r+");
            if (!fichierE) { // if dosen't exist
                cout << "The file " << fichierLecture << " doesn't exist!" << endl; // output dosn't exist
                remove("blocage.txt"); // remove blocking file, we don't need to protect the file that do not exist
            }

            compteur_tentative = 10; // there is no blocking file or we have open it, stop attempt

        } else { // if there is a blocking file
            // wait a little bit
            cout << "Sleep " << compteur_tentative << endl;
            #ifdef __unix__
            usleep(10000);
            #else
            Sleep(10000);
            #endif

            compteur_tentative++;
        }
    }
}

// close parameters file
void fermerFichierE()
{
    fclose(fichierE);

    #ifdef __unix__
    usleep(2000);
    #else
    Sleep(2000);
    #endif

    remove("blocage.txt"); // dont forget to remove the blocking file, we have finish to work with paramaters

}


// Reads parameter values from input file,
// returns 1 if it reaches the end of input file, else returns 0:
bool lireFichier(int &Numr, int &Nr,int &pr, int &max_alleles, double &dpr, double &dgr, double &ar, double &sr, double &hr, double &slr, double &hlr, double &thr, double &gammar,
                     double &Uir, int &itr, double &pasir,double &seuilr, double &Usir, double &Uscr, double &Uner, int &mUscr, int &menager, double &Lr,
                     int &NbPrelimAr, int &NbPrelimBr, int &NbGenr, int &pasprr, int &passcr, int &sampleSr, int &nbr_locus_neutrer, int &position_stream, int &etatr)
{
	int x;
	bool term;
	do {
        x = fgetc(fichierE); //read caracters one by one
	}
	
	
	// until you reach a * or a : or the end of the file
	while ( !( (x == '*') || (x == ':') || (x == EOF) ) ); 
	// * is for a new simulation
	// : is for start from the last inbreeding depression value_comp
	// ! is that we work on this paramaters line
	
	if (x == EOF) {
		cout << "\nEnd of input file\n";
		term = true;
	} else if (x == '*') { // if we should start a new set of simulations

        fseek(fichierE,-1,SEEK_CUR);
        position_stream = ftell(fichierE);
        putc('!',fichierE); // showed that we work on this line in parameters
        fseek(fichierE,0,SEEK_CUR);

        fscanf(fichierE,",%d,",&Numr); // simulation number
        fscanf(fichierE,"%d,",&Nr); // population size
        fscanf(fichierE,"%d,",&pr); // number of deme
        fscanf(fichierE,"%d,",&max_alleles); // maximum number of different S allele
        fscanf(fichierE,"%lf,",&dpr); // pollen dispersal rate
        fscanf(fichierE,"%lf,",&dgr); // seed dispersal rate
        fscanf(fichierE,"%lf,",&ar); // autopolination rate
    	fscanf(fichierE,"%lf,",&sr); // deleterious mutations selection intensity
        fscanf(fichierE,"%lf,",&hr); // deleterious mutations dominance
        fscanf(fichierE,"%lf,",&slr); // lethal mutations selection intensity
        fscanf(fichierE,"%lf,",&hlr); // lethal mutations dominance
        fscanf(fichierE,"%lf,",&thr); // proportion of lethal mutations
        fscanf(fichierE,"%lf,",&gammar); // linkage desequilibrium between S locus and mutations, 0 no link, 1 linkage
		fscanf(fichierE,"%lf,",&Uir); // mutation intensity at the begining
        fscanf(fichierE,"%d,",&itr); // number of step, if >0 search for a threshold if <0 browse U values step by step
        fscanf(fichierE,"%lf,",&pasir); // step width
        fscanf(fichierE,"%lf,",&seuilr); // SC frequency threashold to consider SC invade
        fscanf(fichierE,"%lf,",&Usir); // SI->SI mutation rate
        fscanf(fichierE,"%lf,",&Uscr); // SI->SC mutation rate
        fscanf(fichierE,"%lf,",&Uner); // neutral mutation rate
		fscanf(fichierE,"%d,",&mUscr); // switch of SC mutation: if = 0, only one mutation can occcure at the same time in the metapopulation
        fscanf(fichierE,"%d,",&menager); // switch for delete fixed deleterious or lethal mutations
        fscanf(fichierE,"%lf,",&Lr); // chromosome size / number of crossing over
        fscanf(fichierE,"%d,",&NbPrelimAr); // number of preliminary generation, phase A, only SI->SI mutations
        fscanf(fichierE,"%d,",&NbPrelimBr); // number of preliminary generation, phase B, SI->SI mutations & del/leth mutations
        fscanf(fichierE,"%d,",&NbGenr); // number generation for SC evolution, phase C, SI->SI mutations & del/leth mutations & SI->SC mutations
        fscanf(fichierE,"%d,",&pasprr); // how much generation between writing results in phase A & B
        fscanf(fichierE,"%d,",&passcr);  // how much generation between writing results in phase C
        fscanf(fichierE,"%d,",&sampleSr); // sample size for inbreeding mesuerement 
        fscanf(fichierE,"%d,",&nbr_locus_neutrer); // number of neutral locus 

        etatr = 0; // change this to 0 to say its a new simulation
        term = false;
		
    } else if (x == ':') {

        fseek(fichierE,-1,SEEK_CUR);
        position_stream = ftell(fichierE);
        putc('!',fichierE); // showed that we work on this line in parameters
        fseek(fichierE,0,SEEK_CUR);

        fscanf(fichierE,",%d,",&Numr); // simulation number
        fscanf(fichierE,"%d,",&Nr); // population size
        fscanf(fichierE,"%d,",&pr); // number of deme
        fscanf(fichierE,"%d,",&max_alleles); // maximum number of different S allele
        fscanf(fichierE,"%lf,",&dpr); // pollen dispersal rate
        fscanf(fichierE,"%lf,",&dgr); // seed dispersal rate
        fscanf(fichierE,"%lf,",&ar); // autopolination rate
    	fscanf(fichierE,"%lf,",&sr); // deleterious mutations selection intensity
        fscanf(fichierE,"%lf,",&hr); // deleterious mutations dominance
        fscanf(fichierE,"%lf,",&slr); // lethal mutations selection intensity
        fscanf(fichierE,"%lf,",&hlr); // lethal mutations dominance
        fscanf(fichierE,"%lf,",&thr); // proportion of lethal mutations
        fscanf(fichierE,"%lf,",&gammar); // linkage desequilibrium between S locus and mutations, 0 no link, 1 linkage
		fscanf(fichierE,"%lf,",&Uir); // mutation intensity at the begining
        fscanf(fichierE,"%d,",&itr); // number of step, if >0 search for a threshold if <0 browse U values step by step
        fscanf(fichierE,"%lf,",&pasir); // step width
        fscanf(fichierE,"%lf,",&seuilr); // SC frequency threashold to consider SC invade
        fscanf(fichierE,"%lf,",&Usir); // SI->SI mutations rate
        fscanf(fichierE,"%lf,",&Uscr); // SI->SC mutations rate
        fscanf(fichierE,"%lf,",&Uner); // neutral mutations rate
		fscanf(fichierE,"%d,",&mUscr); // switch of SC mutation: if = 0, only one mutation can occcure at the same time in the metapopulation
        fscanf(fichierE,"%d,",&menager); // switch for delete fixed deleterious or lethals mutations
        fscanf(fichierE,"%lf,",&Lr); // chromosome size / number of crossing over
        fscanf(fichierE,"%d,",&NbPrelimAr); // number of preliminary generation, phase A, only SI->SI mutations
        fscanf(fichierE,"%d,",&NbPrelimBr); // number of preliminary generation, phase B, SI->SI mutations & del/leth mutations
        fscanf(fichierE,"%d,",&NbGenr); // number generation for SC evolution, phase C, SI->SI mutations & del/leth mutations & SI->SC mutations
        fscanf(fichierE,"%d,",&pasprr); // how much generation between writing results in phase A & B
        fscanf(fichierE,"%d,",&passcr);  // how much generation between writing results in phase C
        fscanf(fichierE,"%d,",&sampleSr); // sample size for inbreeding mesuerement 
        fscanf(fichierE,"%d,",&nbr_locus_neutrer); // number of neutral locus 

        etatr = 1; // change this to 1 to say we come back to a started (and stoped for any reason) simulation
        term = false;
    }

	return term;
}
