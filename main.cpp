#include "depression.h"
#include "MersenneTwister.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#ifdef __unix__
#include <sys/stat.h>
#else
#include <direct.h>
#endif
using namespace std;

FILE * fichierP;

// Random number generator:
// MTRand rnd // random seed 
MTRand rnd(1); // choose a seed 

// Pointers on input and output files:
FILE * fichierE;

int main()
{

	// some declaration
    // Parameters:
    int Nt, NbPrelimA, NbPrelimB, NbGen, paspr, passc, sampleS, it, p, mUsc, menage, max_alleles, nbr_locus_neutre;
	double a, s, h, sl, hl, th, gamma, U, pasi, Usi, Usc, Une, L, seuil, Ui, dp, dg;
	result res;
	int position_stream, etat;

	// Opens input and output files:
	bool end;
	end = false;
	bool test;
	bool stop;

	int i;
    int no;
    int numsimul;

	do	{
        //reads parameter values from input file:
        ouvrirFichierE();
 		end = lireFichier(numsimul ,Nt, p, max_alleles, dp, dg, a, s, h, sl, hl, th, gamma,
                    Ui, it, pasi, seuil, Usi, Usc, Une, mUsc, menage, L,
                    NbPrelimA, NbPrelimB, NbGen, paspr, passc, sampleS, nbr_locus_neutre, position_stream, etat);
        fermerFichierE();
		// end = true if end of input file and that there is no more parameters to be processed

        // if population size and deme size do not give a good number of individual by deme, comput a new
        int taille_pop = Nt-Nt%p;

        U = Ui; // number of del/leth mutations start with the initial value (that we protect by a copy)
        stop = false; // initiate stopping system

		if (!end==true) { // if we didn't reached end of paramaters file but we have find paramaters to processed

            if (etat==0) { // if it is a new set of simulations

                // create a directory for output
                char nomdossier[50];
                sprintf(nomdossier,"%d_dossier",numsimul);

                #ifdef __unix__
                mkdir(nomdossier,07777);
                #else
                mkdir(nomdossier);
                #endif

				// create a file just to say that we have begin to work on this simulation and that we have create the directory well
                char nomFichierDebut[256];
                stringstream nomDebut;
                nomDebut << numsimul << "_dossier/" << "debut.txt";
                nomDebut >> nomFichierDebut;

				// if all is good with the directory
                if(fopen(nomFichierDebut,"r")==NULL) {
                    ofstream foutdebut(nomFichierDebut);
                    foutdebut << "Debut" << endl;
                    foutdebut.close();
					
					// print that we processed 
                    cout << "Traitement " << numsimul << endl;

                    // creating a file that summary the differents simulations 
                    char nomFichierParamGeneral[256];
                    stringstream nomFPG;
                    nomFPG << numsimul << "_dossier/" << numsimul << "_paramgeneral.txt";
                    nomFPG >> nomFichierParamGeneral;
                    ofstream foutparamgene(nomFichierParamGeneral);
                    // writing header of column
					// set of simulation number | iteration/simulation number | del/let mutation rate | inbreeding depression at the end | SC frequency at end | result of invasion test | should we stop ? |
                    foutparamgene  << "no" << " " << "iterations" << " " << "U" << " " << "Delta_fin" << " " << "fSC_fin"<< " " << "test" << " " << "stop" << endl;
                    foutparamgene.close();

                    no = 0; // tracking the number of simulation/iteration
                    if (it > 0 ) { // if we look at a threshold value
                        for (i = 0; i < it; i++) { // for the number of simulation/iteration
                            if (i > 0) { // if it is not the first iteration

                                if (res.frSC == 1) { // if SC totaly invade
                                    test=true; // invasion is TRUE
                                    if (res.delta>=0.9900) { // if total invasion and hight inbreeding depression we can stop
                                        stop = true; // stop
                                    }
                                } else if (res.frSC > seuil) { // if SC frequency is higher than SC threeshold in parameters
                                    test = true;  // invasion is TRUE
                                     if (res.delta>=0.995) { // if invasion and very hight inbreeding depression we can stop
                                        stop = true; // stop
                                    }
                                } else { // SC frequency under the expected threshold for an invasion
                                    test = false; // invasion is FALSE
                                    if (res.delta < 0.01) { // if SC dit not invade until there is no inbreeding we can stop
                                        stop = true; // stop
                                    }
                                }

								// write parameters of the previous simulation/iteration
                                foutparamgene.open(nomFichierParamGeneral,std::ofstream::app);
                                foutparamgene << numsimul << " " << no << " " << U << " " << res.delta << " " << res.frSC << " " << test << " " << stop << endl;
                                foutparamgene.close();

                                // if we don't stop, we modifie U (mutation rate)
                                if (stop == false) {
                                    // increase U if no invasion (find a U for which SI maintain)
                                    if ((U == 0) && test) {
                                        U = 1;
                                        i = 0;
                                    } else if ((U == 1) && test) {
                                        U = 2;
                                        i = 0;
                                    } else if ((U == 2) && test) {
                                        U = 3;
                                        i = 0;
                                    } else if ((U == 3) && test) {
                                        U = 4;
                                        i = 0;
                                    } else if ((U == 4) && test) {
                                        U = 5;
                                        i = 0;
                                    } else if ((U == 5) && test) {
                                        U = 6;
                                        i = 0;
                                    } else if ((U == 6) && test) {
                                        U = 7;
                                        i = 0;
                                    } else if ((U == 7) && test) {
                                        U = 8;
                                        i = 0;
                                    } else if ((U == 8) && test) {
                                        U = 9;
                                        i = 0;
                                    } else if ((U == 9) && test) {
                                        U = 10;
                                        i = 0;
                                    }
                                   // if SC did not invade, U = U+( (false-true)/pas^i)
                                   // pas=step is define by parameters 
								   // each new simulation/iteration made the change in U smaller for more and more precision around the U/Inbreeding depression threshold
                                   else ( U += double (test - !test) * pow(pasi,i) );
                                } else {  //stop == true
                                    i = it; // to stop
                                    break; // stop the simulation set
                                }

                            } // end of: if (i > 0) { // if it is not the first iteration

                            // Launch a simulation/iteration
                            cout << no << ", U: " << U << endl;
                            no++;
                            res = recursion(numsimul, no, taille_pop, a, s, h, sl, hl, th, gamma, U, Usi, Usc, Une, mUsc, L, seuil,
                                                NbPrelimA, NbPrelimB, NbGen, paspr, passc, sampleS, nbr_locus_neutre, dp, dg, p, menage, max_alleles);
                        } // end of: for (i = 0; i < it; i++) { // for the number of simulation/iteration
                        cout << endl;
						
						// we have do all the iteration and not because we want to stop, so write the last result
                        if (stop==false) {
                            foutparamgene.open(nomFichierParamGeneral,std::ofstream::app);
                            foutparamgene << numsimul << " " << no << " " << U << " " << res.delta << " " << res.frSC << " " << "NA" << " " << "NA" << endl;
                            foutparamgene.close();
                        }

                    } else if (it < 0) { // if we browse around U values instead of search for a threshold
                        int itp=-it; // change the number of simulation/iteration to a positive number
                        for (i = 0; i < itp; i++) { // for each simulation/iteration
                             no++;
                             // Launch a simulation:
                             cout << no << ", U: " << U << endl;
                             res = recursion(numsimul, no, taille_pop, a, s, h, sl, hl, th, gamma, U, Usi, Usc, Une, mUsc, L, seuil,
                                                NbPrelimA, NbPrelimB, NbGen, paspr, passc, sampleS, nbr_locus_neutre, dp, dg, p, menage, max_alleles);


							// write sumary of results
                            foutparamgene.open(nomFichierParamGeneral,std::ofstream::app);
                            foutparamgene << numsimul << " " << no << " " << U << " " << res.delta << " " << res.frSC << " " << "NA" << " " << "NA" << endl;
                            foutparamgene.close();

                            U += pasi; // increase U for the next one
                        }
                    }

					// we have finish
					// open parameters
                    ouvrirFichierE();
                    fseek(fichierE,position_stream,SEEK_CUR);
                    putc('#',fichierE); // put a # in parameters to say it is finish
                    fermerFichierE();

					// just a file to say it is finish and ok (you can automaticely browse directory and look for this file)
                    char nomFichierFin[256];
                    stringstream nomFinal;
                    nomFinal << numsimul << "_dossier/" << "fin.txt";
                    nomFinal >> nomFichierFin;
                    ofstream foutfin(nomFichierFin);
                    foutfin << "Fini" << endl;
                    foutfin.close();

                    remove(nomFichierDebut); // remove the file to say we have begin

                    cout << endl;
                } else { // there is already a "début"(begining) file
                    cout << "Travail déjà en cours sur " << numsimul << endl; //work in progress
                }
            } else if (etat==1) { // if  etat=1 start again a set of simulation we have already worked on ( ":" at the begining of the row in parameters) 
                cout << "Reprise traitement " << numsimul << endl; // rework

                no = 0; // counter for the number of simulation/iteration

				// declaration for summary file
                char nomFichierParamGeneral[256];
                stringstream nomFPG;
                nomFPG << numsimul << "_dossier/" << numsimul << "_paramgeneral.txt";
                nomFPG >> nomFichierParamGeneral;

                if (it > 0 ) { // we are searching for a threshold
                    cout << "u pivot" << endl; // threshold of mutation rate U

					// open summary
                    fichierP = fopen(nomFichierParamGeneral,"r");

					// read the already initiated summary
                    char testchar[20];
                    for (int k=0;k<7;k++) {
                        fscanf(fichierP,"%s",testchar);
                    }
                    double x;
                    int y=1;
                    int idepart=1;
                    while (y==1) {
                        cout << endl;
                        y = fscanf(fichierP,"%d",&numsimul); // number of the set of simulation 
                        if (y!=-1) { // line by line of the summary file
                            fscanf(fichierP,"%d",&no); // number of simulation/iteration
                            fscanf(fichierP,"%lf",&U); // U value
                            fscanf(fichierP,"%lf",&res.delta); // inbreeding depression
                            fscanf(fichierP,"%lf",&res.frSC); // SC frequency
                            fscanf(fichierP,"%d",&test); // test of invasion
                            fscanf(fichierP,"%d",&stop); // stop or not
							// print what we read
                            cout << numsimul <<  " " << no <<  " " << U << " " << res.delta <<  " " << res.frSC <<  " " << test <<  " " << stop << endl;
                            if (U!=floor(U)) { // if U its not an integer, 
                                idepart++; // it count as an iteration, we tryed to be more precise for the threshold
                            }
                            cout << "idepart: " << idepart << endl;
                        }
                    }

                    int nodepart=no;
                    fclose(fichierP);
					// end of reading and set up, lets start simulation again
                    cout << "fin de lecture, début calcul: " << endl;
                    for (i = idepart; i < it; i++) {
                        
                        if (no > nodepart) {

                            if (res.frSC == 1) { //if last SC frequency was a total invasion
                                test=true; // invasion = TRUE
                                if (res.delta>=0.9900) { // if inbreeding depression was high
                                    stop = true; // stop
                                }
                            } else if (res.frSC > seuil) { // if SC frequency is above the expected threshold for an invasion
                                test = true;  // invasion = TRUE
                                 if (res.delta>=0.9975) { // if inbreeding depression was very high
                                    stop = true; // stop
                                }
                            } else {
                                test = false; // SC frequency was below invasion threshold
                                if (res.delta < 0.01) { // if inbreeding depression is very low
                                    stop = true; // stop
                                }
                            }

							// write the last line of summary
                            fichierP = fopen(nomFichierParamGeneral,"a");
                            fprintf(fichierP, "%i %i %lf %lf %lf %i %i\n",numsimul,no,U,res.delta,res.frSC,test,stop);
                            fclose(fichierP);
                        }


                        // if we dont stop we need to modify U, same process than for a new simulation
                        if (stop == false) {
                            if ((U == 0) && test) {
                                U = 1;
                                i = 0;
                            } else if ((U == 1) && test) {
                                U = 2;
                                i = 0;
                            } else if ((U == 2) && test) {
                                U = 3;
                                i = 0;
                            } else if ((U == 3) && test) {
                                U = 4;
                                i = 0;
                            } else if ((U == 4) && test) {
                                U = 5;
                                i = 0;
                            } else if ((U == 5) && test) {
                                U = 6;
                                i = 0;
                            } else if ((U == 6) && test) {
                                U = 7;
                                i = 0;
                            } else if ((U == 7) && test) {
                                U = 8;
                                i = 0;
                            } else if ((U == 8) && test) {
                                U = 9;
                                i = 0;
                            } else if ((U == 9) && test) {
                                U = 10;
                                i = 0;
                            } else {
								U += double (test - !test) * pow(pasi,i);
							}
                        } else {  //stop == true
                            i = it;
                            break;
                        }

                        // run a simulation/iteration
                        no++;
                        cout <<no << ", U: " << U << endl;
                        res = recursion(numsimul, no, taille_pop, a, s, h, sl, hl, th, gamma, U, Usi, Usc, Une, mUsc, L, seuil,
                                            NbPrelimA, NbPrelimB, NbGen, paspr, passc, sampleS, nbr_locus_neutre, dp, dg, p, menage, max_alleles);
					} // end of: for (i = idepart; i < it; i++) {

					
					if (stop==false) { // if we have done all the simulation without stoping, write the last results in summary
                        fichierP = fopen(nomFichierParamGeneral,"a");
                        fprintf(fichierP, "%i %i %lf %lf %lf NA NA\n",numsimul,no,U,res.delta,res.frSC);
                        fclose(fichierP);
					}

				} else if (it < 0) { // if we browse into U values, same process than before
					cout << "u parcouru" << endl; // browsing
					int itp=-it; // need a positive value

					// open summary for reading
					fichierP = fopen(nomFichierParamGeneral,"r");

					char testchar[20];
					for (int k=0;k<7;k++) {
						fscanf(fichierP,"%s",testchar);
					}
					double x;
					int y=1;
					int idepart=0;
					while (y==1) {
						cout << "\n";
						y = fscanf(fichierP,"%d",&numsimul); // numero de la simul
						if (y!=-1) {
							fscanf(fichierP,"%d",&no); // nbr de valeurs de U testée
							fscanf(fichierP,"%lf",&U); // valeur de U
							fscanf(fichierP,"%lf",&res.delta); // inbreeding
							fscanf(fichierP,"%lf",&res.frSC); // part de SC fin
							fscanf(fichierP,"%s",&test); // test
							fscanf(fichierP,"%s",&stop); // stop
							cout << numsimul <<  " " << no <<  " " << U <<  " " << no <<  " " << res.delta <<  " " << res.frSC <<  " " << test <<  " " << stop << "\n";
							if (U!=floor(U)) {
								idepart++;
							}
						}
					}

					fclose(fichierP); // close summary

					U += pasi;

					// browse from the last U value
					for (i = idepart; i < itp; i++) {
						no++;
                        cout << no << ", U: " << U << endl;
						// run a simulation
                        res = recursion(numsimul, no, taille_pop, a, s, h, sl, hl, th, gamma, U, Usi, Usc, Une, mUsc, L, seuil,
                                        NbPrelimA, NbPrelimB, NbGen, paspr, passc, sampleS, nbr_locus_neutre, dp, dg, p, menage, max_alleles);


						fichierP = fopen(nomFichierParamGeneral,"a");
						fprintf(fichierP, "%i %i %lf %lf %lf NA NA\n",numsimul,no,U,res.delta,res.frSC);
						fclose(fichierP);

						U += pasi;
					}
				}

				// modify parameters to say we are done with this line
				ouvrirFichierE();
				fseek(fichierE,position_stream,SEEK_CUR);
				putc('#',fichierE);
				fermerFichierE();

				// finish file
				char nomFichierFin[256];
				stringstream nomFinal;
				nomFinal << numsimul << "_dossier/" << "fin.txt";
				nomFinal >> nomFichierFin;
				ofstream foutfin(nomFichierFin);
				foutfin << "Fini" << endl;
				foutfin.close();

				// erase start file
				char nomFichierDebut[256];
				stringstream nomDebut;
				nomDebut << numsimul << "_dossier/" << "debut.txt";
				nomDebut >> nomFichierDebut;
				remove(nomFichierDebut);

				cout << endl;

            }
		}

    // end = true if end of input file
	} while (!end);
	return 0 ;
}
