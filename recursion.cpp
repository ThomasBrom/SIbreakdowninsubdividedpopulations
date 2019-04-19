#include "depression.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <cmath>
#include <csignal>
#include <algorithm>
using namespace std;

extern MTRand rnd;
extern FILE * fichierE;


/* function recursion: iterates the life cycle.

   During the first NbPrelimAv generations: Si/Sj mutations, no deleterious mutations, no SC alleles
   During the next NbPrelimBv generations: Si/Sj mutations, deleterious mutations, no SC alleles
   During the next NbGenv generations: Si/Sj mutations + Si/Sc mutations + deleterious mutations

   Other parameters are:
    Numsimul: number of the set of simulation
    iv: number of the simulation/iteration
    Nv: number of diploid individuals
	av: self-pollination rate
	sv, hv: selection and dominance coefficients of deleterious alleles
	slv, hlv: selection and dominance coefficients of lethal alleles
	thv: proportion of lethal mutations
	gammav: location of lethal mutations on the chromosome (0 -> linked to the S-locus, 1 -> over all the chromosome)
	Uv: deleterious mutation rate per haploid genome
	Usiv: mutation rate from Si to Sj
	Uscv: mutation rate from SI to SC and from SC to Si
	Une: mutation rate of enutral allele
	Lv: genome map length (average number of cross-overs at meiosis)
	seuilv: frequency of SC below which SI is considered maintained
	pasprv: number of generations between measures of inbreeding depression and heterosis in preliminary generations
	pasv: number of generations between measures of selfing rates, depression and heterosis
	sampleSv: sample size for estimating inbreeding depression and heterosis
	dpv: pollen dispersal rate
	dgv: seed dispersal rate
	pv: number of deme
*/

result recursion(int numsimulv, int iv, int Nv, double av, double sv, double hv, double slv, double hlv, double thv, double gammav, double Uv,
                     double Usiv, double Uscv, double Unev, int mUscv, double Lv, double seuilv, int NbPrelimAv, int NbPrelimBv, int NbGenv, int pasprv,
                     int pasv, int sampleSv, int nbr_locus_neutrev, double dpv, double dgv, int pv, int menagev, int nbr_allelesv)
{
    //i,j,k,gen: counter, chrmut: mutated chromosome, chrm= maternal chromosome, chrd= paternal chromosome, ind = individual
	unsigned int un_i, un_j, un_k;
	int i, j, k, gen, chrmut, place_chrmut, chrm, chrd, chr1, chr2, chr_foc, chr_comp, chr_glo, chr_loc, ind, nb, nb3, ind_chr, par1, Si, Sj, neutre, mom, dad, SCch, SCchm, SCch_patch, mutSI, mutne, mut, elim_del, elim_let, mut_state;
	int identite_lie_glo, identite_lie_loc, identite_lie_indiv, identite_ind_glo, identite_ind_loc, identite_ind_indiv;

	// some SC frequency and number of fixed mutations erased
	SCch=0 ; SCchm=0 ;SCch_patch=0; elim_let=0; elim_del=0;
	// some fitness values
	double w, w_m, wbar, varw, wmax, wmax_m, le, le_l, rd, fSCbar,slf, depbar, D, D_patch, D_patch_var, div_temp_patch;
	double w_rand, w_out_glo, w_out_glo_SI, w_out_loc, w_out_loc_SI, w_self, reject_glo, reject_loc, reject_fit_glo, reject_fit_loc;
	double w_rand_var, w_out_glo_var, w_out_glo_SI_var, w_out_loc_var, w_out_loc_SI_var, w_self_var;
	result Res;

	// stop signal
	bool signal_d_arret = false;

    chr off1, off2; //chromosome for inbreeding calculation
	vector <int> S_alleles; S_alleles.clear(); //SI allele richness in the population
    vector <int> S_alleles_patch; S_alleles_patch.clear(); // SI allele richness in a deme
    vector <int> SC_alleles; SC_alleles.clear(); // SC allele richness in the population
    vector <int> SC_alleles_patch; SC_alleles_patch.clear(); // SC allele richness in a deme
	vector <double> fSC; fSC.clear(); // SC frequency to determine if we have an invasion
	vector <double> dep; dep.clear(); // Inbreeding depression at the end of the simulation
	vector <double> freq; // frequency of each SI and the SC in metapopulation
	vector <double> freq_loc; // frequency of each SI and the SC in each deme
	freq.resize(nbr_allelesv);
	freq_loc.resize(nbr_allelesv);

	// deme size
	int taille_patch = Nv/pv; // deme size (we have already adjust Nv in "main")
	int two_taille_patch = 2*taille_patch; // number of chromosome in a deme
	int twoN = 2*Nv; // number of chromosome in the metapopulation
	double twoL  = 2*Lv; 

	double Whet = 1 - hv * sv; // fitnes loss for heterozygosity of a deleterious mutation,
	double Whom = 1 - sv; // fitnes loss for homozygosity of a deleterious mutation,
	double Whet_l = 1 - hlv * slv; // fitnes loss for heterozygosity of a lethal mutation,
	double Whom_l = 1 - slv; // fitnes loss for homozygosity of a lethal mutation,

	int cmptSelf = 0; // counter for the number of selfing reproduction
	int cmptgen = 0; // counter for the number of generation with only SC allele (for stop)

	// stuff for deme diversity
	double nS_patch_bar; //SI mean by deme
	double nS_patch_var; //SI variance across deme
    double nSC_patch_bar; //SC mean by deme
	double nSC_patch_var; //SC variance across deme
    double SC_patch_bar; // another SC mean
	double SC_patch_var; // another SC variance
	// dispersal
	double p_disp_bar; // number of pollen dispersal events
    double g_disp_bar; // number of seed dispersal events

	// set up metapopulation
	chr * pop = new chr [twoN]; // metapopulation
   	chr * temp = new chr [twoN]; // temporary metapopulation (for new generation)

	// table of fitness values (for all individuals)
	double * Wij = new double [Nv];
    double * Wij_m = new double [Nv];
    double * fitness_temp = new double [Nv];
	// for each S-allele, sum of fitnesses of compatible chromosomes
	vector<double> Pfec; // for sum of fitnesses
	Pfec.resize(nbr_allelesv+1); // need one for each SI and for SC
    double tab_Pfec[pv][nbr_allelesv+1]; // for sum of fitnesses by deme
    double wbar_p[pv]; // mean fitness by deme

	// creating an outpute file for parameters of the simulation
    char nomFichierParam[256];
	stringstream nomFP;
    nomFP << numsimulv << "_dossier/" << numsimulv << "_param_" << iv <<".txt";
    nomFP >> nomFichierParam;
	ofstream foutparam(nomFichierParam);
    // parameters of the simulation
	foutparam  << "no"<<" " << "i" <<" "<< "N"<< " " << "a"<< " " << "s"<< " " << "h"<< " " << "sl"<< " " << "hl"<< " " << "th"<< " " << "gamma"<< " " << "U"<< " " <<
                     "Usi"<< " " << "nmaxSI"<< " " << "Usc"<< " " << "mUsc"<< " "  << "L"<< " " << "seuil" << " " << "nb_prelim_A"<< " " << "nb_prelim_B"<< " " << "nb_gen"<< " " << "paspr"<< " " <<
                     "pas"<< " " << "nb_sample"<< " " << "nb_locus_neutre"<< " " << "dp"<< " " << "dg"<< " " << "p"<< " " << "menage"<< endl;
    foutparam << numsimulv << " " << iv << " " << Nv << " " << av << " " << sv << " " << hv << " " << slv << " " << hlv << " " << thv << " " << gammav << " " << Uv << " " <<
                     Usiv << " " << nbr_allelesv << " " << Uscv << " " << mUscv << " " << Lv << " " << seuilv << " " << NbPrelimAv << " " << NbPrelimBv << " " << NbGenv << " " << pasprv << " " <<
                     pasv << " " << sampleSv << " " << nbr_locus_neutrev << " " << dpv << " " << dgv << " " << pv << " " << menagev << endl;
    foutparam.close();


    // creating an outpute file for follow simulation every "paspr" then "pas" generation
	char nomFichierResult[256];
	stringstream nomFResult;
	nomFResult << numsimulv << "_dossier/" << numsimulv << "_result_" << iv <<".txt";
	nomFResult >> nomFichierResult;
	ofstream foutresult(nomFichierResult);
	// results
    foutresult <<"generation"
        << " " << "haplotype_number" << " " << "haplotype_diversity"
        << " " << "local_haplotype_number" << " " << "local_haplotype_diversity"  << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance"
        << " " << "SC_frequency" << " " << "selfing_rate" << " " << "SC_mean" << " " << "SC_variance" << " " << "SC_haplotype_number" << " " << "SC_haplotype_diversity" << " " << "SC_haplotype_variance"
        << " " << "fitness_mean" << " " << "fitness_variance" << " " << "mut_per_chr_del" << " " << "mut_per_chr_let" << " " << "fixe_del" << " " << "fix_let"
        << " " << "fit_rand" << " " << "fit_out_glo" << " " << "fit_out_glo_SI" << " " << "fit_out_loc" << " " << "fit_out_loc_SI" << " " << "fit_self"
        << " " << "var_fit_rand" << " " << "var_fit_out_glo" << " " << "var_fit_out_glo_SI" << " " << "var_fit_out_loc" << " " << "var_fit_out_loc_SI" << " " << "var_fit_self"
        << " " << "rejection_glo" << " " << "rejection_loc" << " " << "rejection_fit_glo" << " " << "rejection_fit_loc"
        << " " << "disp_pol" << " " << "disp_graine" <<  endl;
    foutresult.close();

	// creating an outpute file for a more precise following of SC alleles
	char nomFichierSC[256];
	stringstream nomFSC;
	nomFSC << numsimulv << "_dossier/" << numsimulv << "_SC_" << iv <<".txt";
	nomFSC >> nomFichierSC;
	ofstream foutSC(nomFichierSC);
	// generation, SC after mutation, SC in new generation, SC richness, number of deme with SC
    foutSC <<"generation" << " " << "SC_num_mut" << " " << "SC_num_selec" << " " << "SC_rich_mut" << " " << "SC_occupation" << endl;
    foutSC.close();

    // outpute for neutral diversity
	char nomFichierFST[256];
	stringstream nomFFST;
	nomFFST << numsimulv << "_dossier/" << numsimulv << "_FST_" << iv <<".txt";
	nomFFST >> nomFichierFST;
	ofstream foutFST(nomFichierFST);
    // generation, neutral linked in an individual, neutral linked metapopulation (global), neutral linked same deme,... same for independent neutral
	// one for each neutral locus
    foutFST <<"generation";
    for (i=0; i < nbr_locus_neutrev; i++) {
        foutFST << " lie_comp" << i+1 << " lie_glo" << i+1 << " lie_loc" << i+1 << " ind_comp" << i+1 << " ind_glo" << i+1 << " ind_loc" << i+1;
    }
    foutFST << endl;
    foutFST.close();

    // snapshot of the population
	char nomFichierPhoto[256];
	stringstream nomPhoto;
	nomPhoto << numsimulv << "_dossier/" << numsimulv << "_photo_" << iv <<".txt";
	nomPhoto >> nomFichierPhoto;
	ofstream foutPhoto(nomFichierPhoto);
	// phase (A B C), deme, SI (or SC) , fitness
	foutPhoto <<"phase" << " " << "patch" << " " << "haplotype" << " " << "fitness" << endl;
    foutPhoto.close();


    // follow simulation work
	char nomFichiersuivi[256];
	stringstream nomFsuivi;
	nomFsuivi << numsimulv << "_dossier/" << numsimulv << "_suivi_" << iv <<".txt";
	nomFsuivi >> nomFichiersuivi;
	ofstream foutsuivi(nomFichiersuivi);

    //////////// Set up metapopulation
	//// Only heterozygote for SI locus
    ///////////
	foutsuivi << "init" << endl;
	for (i = 0; i < Nv; i++) { // for each individual in the metapopulation
		nb = 2*i; // two chromosomes by individual
		pop[nb].S = rnd.randInt(nbr_allelesv-1) + 1; // random draw of an SI allele
		do { // until we have homozygosity for SI locus
				pop[nb + 1].S = rnd.randInt(nbr_allelesv-1) + 1; // random draw of an SI allele
		} while (pop[nb + 1].S == pop[nb].S); // if its the same than the first, draw again
		// netral set up
		for (j=0;j < nbr_locus_neutrev; j++) {
            pop[nb].neutral_S.push_back(rnd.randInt());
            pop[nb+1].neutral_S.push_back(rnd.randInt());
            pop[nb].neutral_ind.push_back(rnd.randInt());
            pop[nb+1].neutral_ind.push_back(rnd.randInt());
		}
	} // end of set up (for each individual)

	// snapshot after set up
    foutPhoto.open(nomFichierPhoto,std::ofstream::app);
    for (i = 0; i < twoN; i++) {
        foutPhoto << "S" << " " << ((i-(i%two_taille_patch))/two_taille_patch)+1  << " " << pop[i].S << " " << 1 << endl;
    }
    foutPhoto.close();

    ////////////* preliminary generations: no deleterious mutations and no SC alleles
	////////////
    ////////////
	foutsuivi << "phase 1" << endl; // phase A/1
	for (gen = 0; gen < NbPrelimAv; gen++)	{ // for each generation in first phase
		foutsuivi << gen << endl; // print generation in "suivi"
        foutSC.open(nomFichierSC,std::ofstream::app);
		foutSC << gen << " " << 0 << " " << 0 << " " << 0 << endl; // print "SC"
		foutSC.close();
		
		//////////////*
		////mutations
		//////////////

		//// mutations SI->SI
		mutSI = poisdev(twoN * Usiv);  // how many mutations for the entire metapopulation
		for(i = 0; i < mutSI; i++) { // for each mutation
			chrmut = rnd.randInt(twoN - 1); // random draw of the mutated chromosome
			Si = pop[chrmut].S; // keep in mind SI value of the mutant before mutation
			do // until we do not have a different SI allele
            {
                 pop[chrmut].S = rnd.randInt(nbr_allelesv-1) + 1; // random draw of the mutation

            }
			while (pop[chrmut].S == Si); // we need a different one
		} // end of : for each mutation
        //// end of SI -> SI mutations

        //// mutations on independent neutral locus
		// same process than for SI->SI
		mutne = poisdev(twoN * nbr_locus_neutrev * Unev);  // number of mutations
		for(i = 0; i < mutne; i++) { 
			chrmut = rnd.randInt(twoN - 1); 
			place_chrmut = rnd.randInt(nbr_locus_neutrev - 1);
			neutre = pop[chrmut].neutral_ind[place_chrmut];
			do 
            {
                 pop[chrmut].neutral_ind[place_chrmut] = rnd.randInt(); 
            }
			while (pop[chrmut].neutral_ind[place_chrmut] == neutre); 
		} 
      
		//// mutations on linked neutral locus
        // same process than for SI->SI
		mutne = poisdev(twoN * nbr_locus_neutrev * Unev);
		for(i = 0; i < mutne; i++) { 
			chrmut = rnd.randInt(twoN - 1); 
			place_chrmut = rnd.randInt(nbr_locus_neutrev - 1);
			neutre = pop[chrmut].neutral_S[place_chrmut];
			do 
            {
                 pop[chrmut].neutral_S[place_chrmut] = rnd.randInt(); 
            }
			while (pop[chrmut].neutral_S[place_chrmut] == neutre)
		} 


        ////////////// mutation phase end
		//////////////*
		//// mating, selection, new generation
		//////////////
		p_disp_bar=0; // 0 pollen dispersal
		g_disp_bar=0; // 0 seed dispersal 

		//// selection of mothers and fathers
        for (ind = 0; ind < Nv; ind++) { // for each new individual needed 
            ind_chr = 2*ind; 
            // sampling the mother:
            if (rnd.randExc() >= dgv) { // if no seed dispersal, local population, same deme
                chrm = (ind_chr-ind_chr%two_taille_patch)+rnd.randInt(two_taille_patch-1); // random draw of an individual of the same deme
            } else { // if seed dispersal

                do { // until we don't find an individual from another deme
                    chrm = rnd.randInt(twoN-1); // random draw of an individual in the whole metapopulation
                } while (ind_chr-ind_chr%two_taille_patch==chrm-chrm%two_taille_patch); // we need a seed from another deme
            }
			// if seed come from another deme: seed dispersal event (can be include before)
            if (ind_chr-ind_chr%two_taille_patch!=chrm-chrm%two_taille_patch) {
                g_disp_bar++;
            }


			// recombination for mother gamete:
            if (chrm % 2 == 0) // if focal chromosome is on an even position
                rec(temp[ind_chr], pop[chrm], pop[chrm + 1], Lv); // the second chromosome of the mother is on the next odd position
            else // // if focal chromosome is on an odd position
                rec(temp[ind_chr], pop[chrm], pop[chrm - 1], Lv); // the second chromosome of the mother is on the previous even position

            Si = pop[chrm].S; // SI haplotype of mother gamete 
            if(chrm % 2 == 0) // if focal chromosome is on an even position
                Sj = pop[chrm + 1].S; // the second chromosome of the mother is on the next odd position and then the second SI allele is...
            else Sj = pop[chrm - 1].S; // the second chromosome of the mother is on the previous even position and then second SI allele is...


            // sampling the father:
            do { // until we found a compatible pollen
                if (rnd.randExc () >= dpv) { // if it is not a pollen dispersal event
                    chrd = int((chrm-chrm%two_taille_patch)+rnd.randInt(two_taille_patch-1)); //sampling the father into the local deme
                } else { // if it is a dispersal event
                    do { // until it is not a pollen from another deme
                        chrd = rnd.randInt(twoN-1); //sampling the father in the whole metapopulation
                    } while (chrm-chrm%two_taille_patch==chrd-chrd%two_taille_patch); // we need a pollen from another deme
                } // end of pollen dispersal or not
            } while ((pop[chrd].S == Si) || (pop[chrd].S == Sj)); // we need a pollen with a SI allele different of both SI alleles of the mother
            // take note of the pollen dispersal event
            if (chrm-chrm%two_taille_patch!=chrd-chrd%two_taille_patch) {
                p_disp_bar++;
            }

			// recombination for father gamette, same process than for mother
            if (chrd % 2 == 0)
                rec(temp[ind_chr + 1], pop[chrd], pop[chrd + 1], Lv);
            else
                rec(temp[ind_chr + 1], pop[chrd], pop[chrd - 1], Lv);

        } // end of: for (ind = 0; ind < Nv; ind++) { // for each new individual needed 
		
        // count the dispersal rate (number of events/number of individuals)
        p_disp_bar /= Nv;
        g_disp_bar /= Nv;

        //// replace old population by the new (temporary) one
        for (i = 0; i < twoN; i++) // for each chromosome
            pop[i] = temp[i]; // replace pop by temporary pop

        ////////////// end of mating, selection, new pop

		//////////////*
		//// analysis and writing of results
		//////////////

		//// allele S frequency, results, stop the simulation if there is to few S allele
		
		// clear and put at 0
		S_alleles.clear();
		nS_patch_bar=0;
		nS_patch_var=0;
		D_patch=0;
		D_patch_var=0;
		for (un_i = 0; un_i < freq.size(); un_i++)
						freq[un_i] = 0;

		// S allele diversity
		for (i = 0; i < twoN; i++) { // for each chromosome

            Si = pop[i].S; // keep S allele of this chromosome
			
			// if it is the first chromosome of a deme
			// we should clear and put at 0 values for the deme
            if ( (i%two_taille_patch)==0) { 
                S_alleles_patch.clear();
                for (un_i = 0; un_i < freq_loc.size(); un_i++)
                    freq_loc[un_i] = 0;
            }
			// keep local frequency of S allele
            freq_loc[Si-1]++;

			// richness in the deme
            for (un_j = 0; un_j < S_alleles_patch.size(); un_j++) // for each allele already seen in the deme
                if (Si == S_alleles_patch[un_j]) // if the allele have been seen in the patch
                    break; // no need to stay in the loop
                if (un_j == S_alleles_patch.size()) // if we do not found the S allele after have check all the S allele of the deme for the moment
                S_alleles_patch.push_back(Si); // it is a "new" allele for the deme, add it

            // if it is the last chromosome of a deme
			// we should compute some values
            if ( (i%two_taille_patch)==(two_taille_patch-1) ) {  //if it is the last chromosome of a deme
                nS_patch_bar += S_alleles_patch.size(); // for the mean S allele richness
                nS_patch_var += S_alleles_patch.size()*S_alleles_patch.size(); // for the variance in S allele richness

                div_temp_patch = 0; // clear
                for (un_i = 0; un_i < freq_loc.size(); un_i++) { // for each different allele in the deme
                    div_temp_patch += pow((freq_loc[un_i] / two_taille_patch),2); // calculation for diversity
                }
                div_temp_patch = 1/div_temp_patch; // calculation for diversity
                D_patch += div_temp_patch; // for mean diversity
                D_patch_var += div_temp_patch*div_temp_patch; // for variance in diversity (König-Huygens formula)

                if (S_alleles_patch.size() < 3) // if there is less than 3 different SI allele in the deme
                    signal_d_arret=true; // we should stop simulation, no possible mating next generation with two alleles
            } // end of "if it is the last chromosome of a deme"

            // for the whole metapopulation
			// richness of S allele, same process than in a deme
            for (un_j = 0; un_j < S_alleles.size(); un_j++)
                if (Si == S_alleles[un_j])
                    break;
            if (un_j == S_alleles.size())
                S_alleles.push_back(Si);

            freq[Si-1]++; // frequency of each SI allele
        }
        if (S_alleles.size() < 3) // if there is less than 3 different SI allele in the metapopulation (verification in a deme should stop simulation also)
            signal_d_arret=true; // we should stop simulation, no possible mating next generation with two alleles
		
        //// richness by deme
        nS_patch_bar = nS_patch_bar/pv; // for the mean S allele richness
        nS_patch_var = (nS_patch_var/pv)-nS_patch_bar*nS_patch_bar;  // for the variance in S allele richness (König-Huygens)
        // diversity by deme
        D_patch = D_patch/pv; // mean
        D_patch_var = (D_patch_var/pv)-D_patch*D_patch; // variance  (König-Huygens)

        ////diversity at the S-locus or metapopulation
        D = 0;
        for (un_i = 0; un_i < freq.size(); un_i++) // for each possible SI alelle
            D += pow((freq[un_i] / twoN),2); // for simpson index
        D = 1/D; //Gini-Simpson index, probability that two type are different rather than be the same

		// Writes number of S-alleles and effective number of S-alleles in output file every "pasprv" generations:
		if(gen % pasprv == 0) { 
			// some computation on the rejection rate of pollen
            reject_glo = 0;
            reject_loc = 0;
            // average rejection rate in global
			for (j = 0; j < sampleSv; j++) { // for the sample size
				chr1 = rnd.randInt(twoN-1); // choose a focal chromosome (pollen)
				if (chr1 % 2 == 0 ) { // if the focal chr is in even position
                    do {
                        chr2 = rnd.randInt(twoN-1); // sample a random pistil
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ; // be sure that the pistil is not from the same individual than pollen
				} else { // odd position
                    do {
                        chr2 = rnd.randInt(twoN-1); // sample a random pistil
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) ) ;// be sure that the pistil is not from the same individual than pollen
				}
				
				if (chr2 % 2 == 0 ) { // if focal chr of pistil is in even position
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2+1].S) ) { // if pollen have the same S allele than any of the S allele of the pistil individual
                        reject_glo++; // incompatibility
                    }
				} else { // if focal chr of pistil is in odd position
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2-1].S) ) {  // if pollen have the same S allele than any of the S allele of the pistil individual
                        reject_glo++; // incompatibility
                    }
				}
				
			} // end of: for (j = 0; j < sampleSv; j++) { // for the sample size

            // average rejection rate in the same deme, same process than before
			for (j = 0; j < sampleSv; j++) { 
				chr1 = rnd.randInt(twoN-1);
				if (chr1 % 2 == 0 ) {
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1); // sample in the same deme
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ;
				} else { 
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1); // sample in the same deme
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) ) ;
				}
				if (chr2 % 2 == 0 ) {
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2+1].S) ) {
                        reject_loc++;
                    }
				} else {
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2-1].S) ) {
                        reject_loc++;
                    }
				}
			}

            /*foutresult <<"generation"
                << " " << "haplotype_number" << " " << "haplotype_diversity" //ligne alleles totales
                << " " << "local_haplotype_number" << " " << "local_haplotype_diversity"  << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance" //ligne alleles locale
                << " " << "SC_frequency" << " " << "selfing_rate" << " " << "SC_mean" << " " << "SC_variance" << " " << "SC_haplotype_number" << " " << "SC_haplotype_diversity" << " " << "SC_haplotype_variance"  //ligne alleles SC
                << " " << "fitness_mean" << " " << "fitness_variance" << " " << "mut_per_chr_del" << " " << "mut_per_chr_let" << " " << "fixe_del" << " " << "fix_let" // ligne fitness et mutations
                << " " << "fit_rand" << " " << "fit_out_glo" << " " << "fit_out_glo_SI" << " " << "fit_out_loc" << " " << "fit_out_loc_SI" << " " << "fit_self" // ligne fitness théoriques
                << " " << "var_fit_rand" << " " << "var_fit_out_glo" << " " << "var_fit_out_glo_SI" << " " << "var_fit_out_loc" << " " << "var_fit_out_loc_SI" << " " << "var_fit_self" // ligne var fitness théoriques
                << " " << "rejection_glo" << " " << "rejection_loc" << " " << "rejection_fit_glo" << " " << "rejection_fit_loc" // ligne proba rejet
                << " " << "disp_pol" << " " << "disp_graine" <<  endl; // ligne disp
            foutresult.close();*/

            foutresult.open(nomFichierResult,std::ofstream::app);
			foutresult << gen
                << " " << S_alleles.size() << " " << D // metapopulation SI allele line
                << " " << nS_patch_bar << " " << D_patch << " " << nS_patch_var << " " << D_patch_var // locale allele line
                << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0  // SC allele line
                << " " << 1 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 // fitness & mutation line
                << " " << 1 << " " << 1 << " " << 1 << " " << 1 << " " << 1 << " " << 1 // theoretical fitness line
                << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 // variance for theoretical fitness
                << " " << reject_glo/sampleSv << " " << reject_loc/sampleSv << " " << "NA" << " " << "NA" // incompatibility probability line
                << " " << p_disp_bar << " " << g_disp_bar << endl; // dispersal line
            foutresult.close();

            // neutral diversity/identity
			foutFST.open(nomFichierFST,std::ofstream::app);
            foutFST << gen;
            for (i = 0; i<nbr_locus_neutrev; i++) { // for each neutral locus
                identite_lie_glo = 0;
                identite_lie_loc = 0;
                identite_lie_indiv = 0;
                identite_ind_glo = 0;
                identite_ind_loc = 0;
                identite_ind_indiv = 0;
                for (j = 0; j < sampleSv; j++) { // for sample size
                    chr_foc = rnd.randInt(twoN-1); // sample a focal chromosome
                    if (chr_foc % 2 == 0 ) { // even position
						// complementary chromosome is:
                        chr_comp = chr_foc+1;
						
                        do { // until it is another individual than the focal one
                            chr_glo = rnd.randInt(twoN-1); // sample an individual in the metapopulation
                        } while ( (chr_glo==chr_foc) | (chr_glo==chr_comp) ) ; // until it is another individual than the focal one
						
                        do { // until it is another individual than the focal one
                            chr_loc = (chr_foc-chr_foc%two_taille_patch)+rnd.randInt(two_taille_patch-1); // sample an individual in the same deme
                        } while ( (chr_loc==chr_foc) | (chr_loc==chr_comp) ) ; // until it is another individual than the focal one

                    } else { // odd position
						// complementary chromosome is:
                        chr_comp = chr_foc-1;
						
                        do { // until it is another individual than the focal one
                            chr_glo = rnd.randInt(twoN-1); // sample an individual in the metapopulation
                        } while ( (chr_glo==chr_foc) | (chr_glo==chr_comp) ) ; // until it is another individual than the focal one
						
                        do { // until it is another individual than the focal one
                            chr_loc = (chr_foc-chr_foc%two_taille_patch)+rnd.randInt(two_taille_patch-1); // sample an individual in the same deme
                        } while ( (chr_loc==chr_foc) | (chr_loc==chr_comp) ) ; // until it is another individual than the focal one

                    }
					// if focal chr have the same linked neutral allele than complementary chr
                    if ( pop[chr_foc].neutral_S[i]==pop[chr_comp].neutral_S[i] ) {
                        identite_lie_indiv++;
                    }
					// if focal chr have the same linked neutral allele than random chr of the metapopulation
                    if ( pop[chr_foc].neutral_S[i]==pop[chr_glo].neutral_S[i] ) {
                        identite_lie_glo++;
                    }
					// if focal chr have the same linked neutral allele than random chr of the same deme
                    if ( pop[chr_foc].neutral_S[i]==pop[chr_loc].neutral_S[i] ) {
                        identite_lie_loc++;
                    }
					// if focal chr have the same independent neutral allele than complementary chr
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_comp].neutral_ind[i] ) {
                        identite_ind_indiv++;
                    }
					// if focal chr have the same independent neutral allele than random chr of the metapopulation
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_glo].neutral_ind[i] ) {
                        identite_ind_glo++;
                    }
					// if focal chr have the same independent neutral allele than random chr of the same deme
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_loc].neutral_ind[i] ) {
                        identite_ind_loc++;
                    }
                } // end of: for (j = 0; j < sampleSv; j++) { // for sample size
			// write results
            foutFST << " " << identite_lie_indiv << " " << identite_lie_glo << " " << identite_lie_loc << " " << identite_ind_indiv << " " << identite_ind_glo << " " << identite_ind_loc;
            } // end of: for each neutral locus
            foutFST << endl;
			foutFST.close();

		} // end of writing results module (every xxx generation)

		////////////// end of analysis and results module

	} // end of phase 1/A

	// snapshot after phase A
    foutPhoto.open(nomFichierPhoto,std::ofstream::app);
    for (i = 0; i < twoN; i++) {
        foutPhoto << "A" << " " << ((i-(i%two_taille_patch))/two_taille_patch)+1 << " " << pop[i].S << " " << 1 << endl;
    }
    foutPhoto.close();

    ////////////
	//////////// end of phase A/1
    ////////////
    if (signal_d_arret==false) // if we can continue

    //////////// 
	//////////// start phase B/2 introduce deleterious/lethals mutations and thus differences of fitness
    ////////////
    foutsuivi << endl;
	foutsuivi <<  "phase 2" << endl;
	for (gen = 0; gen < NbPrelimBv; gen++) { // for each generation of the phase B

        foutsuivi << gen << endl;
        foutSC.open(nomFichierSC,std::ofstream::app);
        foutSC << gen+NbPrelimAv << " " << 0 << " " << 0 << " " << 0 << endl; // still no SC allele
        foutSC.close();
        //////////////**
		////mutations
		//////////////

		//// mutations SI->SI, same than in phase A
        mutSI = poisdev(twoN * Usiv);

		for(i = 0; i < mutSI; i++) { 
			chrmut = rnd.randInt(twoN - 1); 
			Si = pop[chrmut].S; 
			do {
                pop[chrmut].S = rnd.randInt(nbr_allelesv-1) + 1;
            }
			while (pop[chrmut].S == Si);
		}

        //// neutral mutations, same as before in phase A
		
		mutne = poisdev(twoN * nbr_locus_neutrev * Unev);  
		for(i = 0; i < mutne; i++) {
			chrmut = rnd.randInt(twoN - 1); 
			place_chrmut = rnd.randInt(nbr_locus_neutrev - 1);
			neutre = pop[chrmut].neutral_ind[place_chrmut]; 
			do 
            {
                 pop[chrmut].neutral_ind[place_chrmut] = rnd.randInt(); 
            }
			while (pop[chrmut].neutral_ind[place_chrmut] == neutre); 
		} 
        
      	mutne = poisdev(twoN * nbr_locus_neutrev * Unev);
		for(i = 0; i < mutne; i++) {
			chrmut = rnd.randInt(twoN - 1);
			place_chrmut = rnd.randInt(nbr_locus_neutrev - 1);
			neutre = pop[chrmut].neutral_S[place_chrmut];
			do
            {
                 pop[chrmut].neutral_S[place_chrmut] = rnd.randInt();
            }
			while (pop[chrmut].neutral_S[place_chrmut] == neutre);
		}

		//// lethal and deleterious mutations:
		for (i = 0; i < twoN; i++) { // for each chromosome

            mut = poisdev(Uv * thv); // number of new lethal mutations on chromosome Uv=mutation frequency, thv=lethal mutation proportion
            for (j = 0; j < mut; j++) { // for each lethal mutation
                //each mutation has a random position between (1-g)L and (1+g)L
                pop[i].let.push_back(Lv * (1 - gammav + (2 * gammav * rnd.rand())));
                sort(pop[i].let.begin(), pop[i].let.end()); // order lethal mutations
            } 

            mut = poisdev(Uv * (1 - thv)); // number of new deleterious mutations Uv=mutation frequency, thv=lethal mutation proportion
            for (j = 0; j < mut; j++) { //for each deleterious mutation
                //each mutation has a random position between 0 and 2L
                pop[i].sel.push_back(twoL * rnd.rand()); 
                sort(pop[i].sel.begin(), pop[i].sel.end()); // order deleterious mutations
            }
        }
        //// end of mutation phase

		//////////////**
		//// selection/mating -> population of the next generation
		//////////////
		
		//// we need to know fitness of each individual
		wbar = 0; // mean fitness
		wmax = 0; // maximum fitness
		varw = 0; // fitness variance
        for (j = 0; j < Nv; j++) { // for each individual
            nb = 2*j; // "first" chromosome, in even position, of the individual is nb
            w = fitness(pop[nb], pop[nb + 1], Whet, Whom, Whet_l, Whom_l); // fitness calculation
            Wij[j] = w; // keep fitness
            wbar += w; // for mean fitness
            varw += w * w; // for variance of fitness
            if (wmax < w) // if this individual have the highest fitness (for the moment)
                wmax = w; // it is the new highest fitness
        }
		varw /= Nv; // for fitness variance
		wbar /= Nv; // for mean fitness
		varw -= wbar * wbar; // for fitness variance (König-Huygens)

		// for keep dispersal event
		p_disp_bar = 0;
        g_disp_bar = 0;

		//// next generation, same process than before, but we now have difference in fitness
        for (ind = 0; ind < Nv; ind++) { // for each individual at the next generation
            ind_chr = 2*ind; 
            // sampling the mother:
            do{ // sample undividuals with the highest fitness
                if (rnd.randExc () >= dgv) { // no seed dispersal
                    chrm = (ind_chr-ind_chr%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                    mom = chrm/2; 

                } else { // seed dispersal

                    do{ 
                        chrm = rnd.randInt(twoN-1); 
                        mom = chrm/2; // individual number (for fitness value)
                    } while (ind_chr-ind_chr%two_taille_patch==chrm-chrm%two_taille_patch);
                }

            } while (rnd.rand() >= Wij[mom] / wmax); // fitness selection ! 
            
			// if seed dispersal
			if (ind_chr-ind_chr%two_taille_patch!=chrm-chrm%two_taille_patch) {
                g_disp_bar++;
            }

			// recombination for the mother gamete, as before:
            if (chrm % 2 == 0)
                rec(temp[ind_chr], pop[chrm], pop[chrm + 1], Lv);
            else 
                rec(temp[ind_chr], pop[chrm], pop[chrm - 1], Lv);

            // maternal S-alleles
            Si = pop[chrm].S; //first S allele
            if(chrm % 2 == 0) // if first S allele/chromosome is in even or odd position:
                Sj = pop[chrm+1].S; // this one is the second
            else Sj = pop[chrm-1].S; // else this one is the second

			
            // sampling the father:
            do { // we need a pollen compatible with mother (Si Sj) and we select on fitness
                if (rnd.randExc () >= dpv) { // if no pollen dispersal
                    chrd = (chrm-chrm%two_taille_patch)+rnd.randInt(two_taille_patch-1); 
                    dad = chrd/2; // to get fitness later
                } else { // if pollen dispersal
                    do { // until we get a pollen from another deme
                        chrd = rnd.randInt(twoN-1); 
                        dad = chrd/2; // to get fitness later
                    } while (chrm-chrm%two_taille_patch==chrd-chrd%two_taille_patch); // need a different deme
                }
			// not the same individual, selection on fitness, compatibility
            } while ((mom == dad)||(rnd.rand() >= (Wij[dad] / wmax) )||(pop[chrd].S == Si)||(pop[chrd].S == Sj)); 
						
            // if pollen dispersal..
            if (chrm-chrm%two_taille_patch!=chrd-chrd%two_taille_patch) {
                p_disp_bar++;
            }
			
            // recombination for the father gamete, as before:
            if (chrd % 2 == 0)
                rec(temp[ind_chr + 1], pop[chrd], pop[chrd + 1], Lv);
            else
                rec(temp[ind_chr + 1], pop[chrd], pop[chrd - 1], Lv);
        } // end of: for (ind = 0; ind < Nv; ind++) { // for each individual at the next generation

        // dispersal rates
        p_disp_bar /= Nv;
        g_disp_bar /= Nv;

		// new population, population of the next generation take place
		for (i = 0; i < twoN; i++)
			pop[i] = temp[i];

        ////////////// end of population creation
				
		//////////////
		//// analysis and writing results
		//////////////
		
		//// same process than in phase A
		S_alleles.clear(); 
		nS_patch_bar=0;
		nS_patch_var=0;
		D_patch=0;
		D_patch_var=0;
		for (un_i = 0; un_i < freq.size(); un_i++)
						freq[un_i] = 0; 

		for (i = 0; i < twoN; i++) { // for each chromosome

            Si = pop[i].S;
            if ( (i%two_taille_patch)==0) { // if first chromosome of the deme
                S_alleles_patch.clear(); 
                for (un_i = 0; un_i < freq_loc.size(); un_i++)
                    freq_loc[un_i] = 0; 
            }
            freq_loc[Si-1]++;

			// S allele richness in the deme
            for (un_j = 0; un_j < S_alleles_patch.size(); un_j++)
                if (Si == S_alleles_patch[un_j])
                    break;
                if (un_j == S_alleles_patch.size()) 
                S_alleles_patch.push_back(Si);

            
            if ( (i%two_taille_patch)==(two_taille_patch-1) ) { // if last chromosome of the deme
                nS_patch_bar += S_alleles_patch.size(); 
                nS_patch_var += S_alleles_patch.size()*S_alleles_patch.size();

                div_temp_patch = 0;
                for (un_i = 0; un_i < freq_loc.size(); un_i++) { 
                    div_temp_patch += pow((freq_loc[un_i] / two_taille_patch),2);
                }
                div_temp_patch = 1/div_temp_patch; 
                D_patch += div_temp_patch; 
                D_patch_var += div_temp_patch*div_temp_patch; 

                if (S_alleles_patch.size() < 3) 
                    signal_d_arret=true;
            } 

            // SI allele richness for metapopulation, as in phase A
            for (un_j = 0; un_j < S_alleles.size(); un_j++)
                if (Si == S_alleles[un_j])
                    break; // 
            if (un_j == S_alleles.size()) 
                S_alleles.push_back(Si);
			// SI allele frequency
            freq[Si-1]++;
			
        } // end of: for (i = 0; i < twoN; i++) { // for each chromosome
		
		// stop simulation if less than 3 different SI alleles
        if (S_alleles.size() < 3) 
            signal_d_arret=true;
		
        
		// some date on SI allele as in phase A
        nS_patch_bar = nS_patch_bar/pv; 
        nS_patch_var = (nS_patch_var/pv)-nS_patch_bar*nS_patch_bar;
        D_patch = D_patch/pv;
        D_patch_var = (D_patch_var/pv)-D_patch*D_patch;

        //// diversity at the S-locus, as before
        D = 0;
        for (un_i = 0; un_i < freq.size(); un_i++)
            D += pow((freq[un_i] / twoN),2); 
        D = 1/D;

		// if we should write results this generation
		if (gen % pasprv == 0) { 
			// number of mutations per chromosome:
			le = 0; // for deleterious mutations
			for (i = 0; i < twoN; i++) // for each chromosome
				le += temp[i].sel.size(); // add the number of mutation

			le_l = 0; // for lethal mutations
			for (i = 0; i < twoN; i++) // for each chromosome
				le_l += temp[i].let.size(); // add the number of mutation

			// if we should erase fixed mutations,
            if (menagev==1) {
                elim_let = menage_let(temp,twoN); // function for lethal mutations
                elim_del = menage_del(temp,twoN); // function for deleterious mutations
            }

			// measuring the average fitness individuals from random mating, outcrossing, selfing....
			w_rand = 0; //fitness for an individual from a totaly random mating
			w_out_glo = 0; // fitness for outcrossing in the whole metapopulation
			w_out_glo_SI = 0; // fitness for outcrossing in the whole metapopulation taking into account SI
			w_out_loc = 0; // fitness for outcrossing in the same deme
			w_out_loc_SI = 0; // fitness for outcrossing in the same deme taking into account SI
			w_self = 0; // fitness for individuals from selfing

			// random fitness
			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1); // mother
				chr2 = rnd.randInt(twoN-1); // father
                if (chr1 % 2 == 0 ) {
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv); 
				} else { 
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv); /
				}
				if (chr2 % 2 == 0) { 
					rec(off2, pop[chr2], pop[chr2 + 1], Lv); 
				} else {
					rec(off2, pop[chr2], pop[chr2 - 1], Lv);
                }
				w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l); // fitness
				w_rand += w_m;
				w_rand_var += w_m*w_m;
			}
			w_rand = w_rand/sampleSv;
			w_rand_var = (w_rand_var/sampleSv)-(w_rand*w_rand);

			// fitness for outcrossing in the metapopulation whithout SI
			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1); // mother
                if (chr1 % 2 == 0 ) { 
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv);
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ;
				} else {
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) ) ;
				}
				if (chr2 % 2 == 0) 
					rec(off2, pop[chr2], pop[chr2 + 1], Lv);
				else
					rec(off2, pop[chr2], pop[chr2 - 1], Lv);
                w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l); // fitness
				w_out_glo += w_m;
				w_out_glo_var += w_m*w_m;
			}
			w_out_glo = w_out_glo/sampleSv;
			w_out_glo_var = (w_out_glo_var/sampleSv)-(w_out_glo*w_out_glo);

            // fitness for outcrossing in the metapopulation whith SI
			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1); // mother
                if (chr1 % 2 == 0 ) { 
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv); 
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1+1==chr2) | (pop[chr1].S == pop[chr2].S) | (pop[chr1+1].S == pop[chr2].S) ) ;
				} else {
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1-1==chr2) | (pop[chr1].S == pop[chr2].S) | (pop[chr1-1].S == pop[chr2].S) ) ;
				}
				if (chr2 % 2 == 0) 
					rec(off2, pop[chr2], pop[chr2 + 1], Lv);
				else
					rec(off2, pop[chr2], pop[chr2 - 1], Lv);
                w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l);
				w_out_glo_SI += w_m;
				w_out_glo_SI_var += w_m*w_m;
			}
			w_out_glo_SI = w_out_glo_SI/sampleSv;
			w_out_glo_SI_var = (w_out_glo_SI_var/sampleSv)-(w_out_glo_SI*w_out_glo_SI);

			// fitness for local outcrossing (same deme) without SI
			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1); // mother
                if (chr1 % 2 == 0 ) { 
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv);
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ;
				} else {
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) ) ;
				}
				if (chr2 % 2 == 0)
					rec(off2, pop[chr2], pop[chr2 + 1], Lv);
				else
					rec(off2, pop[chr2], pop[chr2 - 1], Lv);
				w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l); // fitness
				w_out_loc += w_m;
				w_out_loc_var += w_m*w_m;
			}
			w_out_loc = w_out_loc/sampleSv;
			w_out_loc_var = (w_out_loc_var/sampleSv)-(w_out_loc*w_out_loc);


			// fitness for local outcrossing (same deme) whith SI
			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1); //tirage de la mere
                if (chr1 % 2 == 0 ) { //si on est en position paire dans le vecteur
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv); //recombinaison pour la mere
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1); //tirage du pere
                    } while ( (chr1==chr2) | (chr1+1==chr2) | (pop[chr1].S == pop[chr2].S) | (pop[chr1+1].S == pop[chr2].S) ) ;
				} else { //si on est en position impaire dans le vecteur
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv); //recombinaison pour la mere
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1); //tirage du pere
                    } while ( (chr1==chr2) | (chr1-1==chr2) | (pop[chr1].S == pop[chr2].S) | (pop[chr1-1].S == pop[chr2].S) ) ;
				}
				if (chr2 % 2 == 0) //on prend le 2ème chromosome de l'individu (+1 ou -1 selon la position dans le vecteur)
					rec(off2, pop[chr2], pop[chr2 + 1], Lv); //recombinaison pour le pere
				else
					rec(off2, pop[chr2], pop[chr2 - 1], Lv); //recombinaison pour le pere
                w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l); //calcule de fitness selfing et ajout
				w_out_loc_SI += w_m;
				w_out_loc_SI_var += w_m*w_m;
			}
			w_out_loc_SI = w_out_loc_SI/sampleSv;
			w_out_loc_SI_var = (w_out_loc_SI_var/sampleSv)-(w_out_loc_SI*w_out_loc_SI);


            // fitness for individuals from selfing
			for (j = 0; j < sampleSv; j++) { 
				chr1 = rnd.randInt(twoN-1); 
				if (chr1 % 2 == 0) { 
					nb3 = chr1 + rnd.randInt(1); 
					rec(off1, pop[chr1], pop[chr1 + 1], Lv); 
				} else { 
					nb3 = chr1 - rnd.randInt(1);
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
				}
				if (nb3 % 2 == 0) 
					rec(off2, pop[nb3], pop[nb3 + 1], Lv);
				else
					rec(off2, pop[nb3], pop[nb3 - 1], Lv); 
				w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l); // fitness
				w_self += w_m;
				w_self_var += w_m*w_m;
			}
			w_self = w_self/sampleSv;
			w_self_var = (w_self_var/sampleSv)-(w_self*w_self);

            // some computation on the rejection rate of pollen
			// we add rejection rate taking into account fitness, otherwise it is the same process
            reject_glo = 0; // metapopulation 
            reject_fit_glo = 0; // metapopulation with fitness
            reject_loc = 0; // in the same deme
            reject_fit_loc = 0; // in the same deme with fitness

			// we should compute fitnesses of the new population
            wmax = 0;
            for (j = 0; j < Nv; j++) {
                nb = 2*j; 
                w_m = fitness(pop[nb], pop[nb + 1], Whet, Whom, Whet_l, Whom_l); 
                Wij_m[j] = w_m;
            if (wmax_m < w_m) 
                wmax_m = w_m; 
            }

            // average rejection rate in global
			for (j = 0; j < sampleSv; j++) { 
				chr1 = rnd.randInt(twoN-1); 
				if (chr1 % 2 == 0 ) { 
                    do {
                        chr2 = rnd.randInt(twoN-1); 
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ;
				} else {
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) ) ;
				}
				if (chr2 % 2 == 0 ) { 
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2+1].S) ) {
                        reject_glo++;
                    }
				} else { 
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2-1].S) ) {
                        reject_glo++;
                    }
				}
			}
            // average rejection rate in global taking into account fitness
			for (j = 0; j < sampleSv; j++) { 
                do {
                    chr1 = rnd.randInt(twoN-1); 
                } while (rnd.rand() >= Wij_m[chr1/2] / wmax_m);
				if (chr1 % 2 == 0 ) {
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1+1==chr2) | (rnd.rand() >= Wij_m[chr2/2] / wmax_m) ) ;
				} else { 
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) | (rnd.rand() >= Wij_m[chr2/2] / wmax_m) ) ;
				}
				if (chr2 % 2 == 0 ) { 
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2+1].S) ) {
                        reject_fit_glo++;
                    }
				} else {
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2-1].S) ) {
                        reject_fit_glo++;
                    }
				}
			}

            // average rejection rate in local/same deme
			for (j = 0; j < sampleSv; j++) { 
				chr1 = rnd.randInt(twoN-1); 
				if (chr1 % 2 == 0 ) {
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1); 
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ;
				} else { 
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1); 
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) ) ;
				}
				if (chr2 % 2 == 0 ) {
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2+1].S) ) {
                        reject_loc++;
                    }
				} else { 
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2-1].S) ) {
                        reject_loc++;
                    }
				}
			}

            // average rejection rate in local taking into account fitness
			for (j = 0; j < sampleSv; j++) {
                do {
                    chr1 = rnd.randInt(twoN-1);
                } while (rnd.rand() >= Wij_m[chr1/2] / wmax_m);
				if (chr1 % 2 == 0 ) { 
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1+1==chr2) | (rnd.rand() >= Wij_m[chr2/2]/wmax_m) ) ;
				} else { 
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1-1==chr2) | (rnd.rand() >= Wij_m[chr2/2]/wmax_m) ) ;
				}
				if (chr2 % 2 == 0 ) { 
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2+1].S) ) {
                        reject_fit_loc++;
                    }
				} else {
                    if ( (pop[chr1].S==pop[chr2].S) || (pop[chr1].S==pop[chr2-1].S) ) {
                        reject_fit_loc++;
                    }
				}
			}

           
            /*foutresult <<"generation"
                << " " << "haplotype_number" << " " << "haplotype_diversity" //ligne alleles totales
                << " " << "local_haplotype_number" << " " << "local_haplotype_diversity"  << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance" //ligne alleles locale
                << " " << "SC_frequency" << " " << "selfing_rate" << " " << "SC_mean" << " " << "SC_variance" << " " << "SC_haplotype_number" << " " << "SC_haplotype_diversity" << " " << "SC_haplotype_variance"  //ligne alleles SC
                << " " << "fitness_mean" << " " << "fitness_variance" << " " << "mut_per_chr_del" << " " << "mut_per_chr_let" << " " << "fixe_del" << " " << "fix_let" // ligne fitness et mutations
                << " " << "fit_rand" << " " << "fit_out_glo" << " " << "fit_out_glo_SI" << " " << "fit_out_loc" << " " << "fit_out_loc_SI" << " " << "fit_self" // ligne fitness théoriques
                << " " << "var_fit_rand" << " " << "var_fit_out_glo" << " " << "var_fit_out_glo_SI" << " " << "var_fit_out_loc" << " " << "var_fit_out_loc_SI" << " " << "var_fit_self" // ligne var fitness théoriques
                << " " << "rejection_glo" << " " << "rejection_loc" << " " << "rejection_fit_glo" << " " << "rejection_fit_loc" // ligne proba rejet
                << " " << "disp_pol" << " " << "disp_graine" <<  endl; // ligne disp
            foutresult.close();*/

			
			// writing result
            foutresult.open(nomFichierResult,std::ofstream::app);
			foutresult << NbPrelimAv+gen
                << " " <<S_alleles.size() << " " << D //ligne alleles totales
                << " " << nS_patch_bar << " " << D_patch << " " << nS_patch_var << " " << D_patch_var //ligne alleles locale
                << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0  //ligne alleles SC
                << " " << wbar << " " << varw << " " << le / twoN<< " " << le_l / twoN << " " << elim_del << " " << elim_let // ligne fitness et mutations
                << " " << w_rand << " "  << w_out_glo << " "  << w_out_glo_SI << " "  << w_out_loc << " "  << w_out_loc_SI << " "  << w_self // ligne fitness théoriques
                << " " << w_rand_var << " "  << w_out_glo_var << " "  << w_out_glo_SI_var << " "  << w_out_loc_var << " "  << w_out_loc_SI_var << " "  << w_self_var // ligne var fitness théoriques
                << " " << reject_glo/sampleSv << " " << reject_loc/sampleSv << " " << reject_fit_glo/sampleSv << " " << reject_fit_loc/sampleSv // ligne proba rejet
                << " " << p_disp_bar << " " << g_disp_bar << endl; // ligne disp
            foutresult.close();

			// write inbreeding depression
		    if (gen >= pasv) { 
                dep.push_back(1 - (w_self / w_out_loc)); 
				//dep.push_back(1 - (wself / wout)); 
		    }

            // neutrals alleles, same than in phase A
			foutFST.open(nomFichierFST,std::ofstream::app);
            foutFST << gen;
            for (i = 0; i<nbr_locus_neutrev; i++) {
                identite_lie_glo = 0;
                identite_lie_loc = 0;
                identite_lie_indiv = 0;
                identite_ind_glo = 0;
                identite_ind_loc = 0;
                identite_ind_indiv = 0;
                for (j = 0; j < sampleSv; j++) {
                    chr_foc = rnd.randInt(twoN-1);
                    if (chr_foc % 2 == 0 ) { 

                        chr_comp = chr_foc+1;
                        do { 
                            chr_glo = rnd.randInt(twoN-1);
                        } while ( (chr_glo==chr_foc) | (chr_glo==chr_comp) ) ;
                        do { 
                            chr_loc = (chr_foc-chr_foc%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                        } while ( (chr_loc==chr_foc) | (chr_loc==chr_comp) ) ;

                    } else { 
                      
                        chr_comp = chr_foc-1;
                        do { 
                            chr_glo = rnd.randInt(twoN-1);
                        } while ( (chr_glo==chr_foc) | (chr_glo==chr_comp) ) ;
                        do {
                            chr_loc = (chr_foc-chr_foc%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                        } while ( (chr_loc==chr_foc) | (chr_loc==chr_comp) ) ;

                    }

                    if ( pop[chr_foc].neutral_S[i]==pop[chr_comp].neutral_S[i] ) {
                        identite_lie_indiv++;
                    }
                    if ( pop[chr_foc].neutral_S[i]==pop[chr_glo].neutral_S[i] ) {
                        identite_lie_glo++;
                    }
                    if ( pop[chr_foc].neutral_S[i]==pop[chr_loc].neutral_S[i] ) {
                        identite_lie_loc++;
                    }
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_comp].neutral_ind[i] ) {
                        identite_ind_indiv++;
                    }
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_glo].neutral_ind[i] ) {
                        identite_ind_glo++;
                    }
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_loc].neutral_ind[i] ) {
                        identite_ind_loc++;
                    }
                } 
            foutFST << " " << identite_lie_indiv << " " << identite_lie_glo << " " << identite_lie_loc << " " << identite_ind_indiv << " " << identite_ind_glo << " " << identite_ind_loc;
            } 
            foutFST << endl;
			foutFST.close();


		} // end of:	if (gen % pasprv == 0) { // if we should write results this generation
		

		// should we stop (less than 3 SI alleles)
		if (signal_d_arret==true) { 
			foutsuivi << "arret ";
			break;
		}

	}
    ////////////
	//////////// 
    ////////////

	// compute fitness for snapshoot
    for (i = 0; i < Nv; i++) {
         fitness_temp[i] = fitness(pop[i], pop[i+1], Whet, Whom, Whet_l, Whom_l);
    }

	// snapshot for phase B
    foutPhoto.open(nomFichierPhoto,std::ofstream::app);
    for (i = 0; i < twoN; i++) {
        foutPhoto << "B" << " " << ((i-(i%two_taille_patch))/two_taille_patch)+1  << " " << pop[i].S << " " << fitness_temp[i/2] << endl;
    }
    foutPhoto.close();


    if (!signal_d_arret==true) // if we can continue

    ////////////
	//////////// last phase (C) introduce SC allele and evolution of the mating system
    ////////////
    foutsuivi << endl;
	foutsuivi << "phase 3" << endl;
	
    for ( gen = 0; gen<NbGenv; gen++) { // for each generation until the end

        foutsuivi << gen << endl;
        foutSC.open(nomFichierSC,std::ofstream::app);
        foutSC << NbPrelimAv+NbPrelimBv+gen << " ";
        foutSC.close();

        //////////////
		////mutations, same process than before but SI->SC mutations
		//////////////
		SCchm = 0; // for keep the number of SC mutations
		mut_state=0; // for know if there is already a mutation in this generation
        for (i = 0; i < twoN; i++) { // for each chromosome

            //// mutation on locus S
            rd = rnd.rand(); // random [0,1]
            if (pop[i].S > 0) { // if it is a SI allele

				if (rd < Uscv) { // if mutation SI->SC
					// if multiple different mutations at the same time is true (=1) or there is no SC allele previously and still no mutation at this generation
                    if (mUscv==1 || (SCch==0 & mut_state==0)) { 
						pop[i].S = -1*pop[i].S; // there is a mutation
						mut_state = 1; // a SC mutation at this generation
						SCchm++; // on more SC mutation
                    }

				} else if (rd < (Uscv + Usiv)) { //if SI->SI
					// same process
					Si = pop[i].S; 
					do pop[i].S = rnd.randInt(nbr_allelesv-1) + 1; 
					while (pop[i].S == Si); 
                }
			}

            //// lethals, same process
            mut = poisdev(Uv * thv); // number of new lethal mutations on chromosome
            for (j = 0; j < mut; j++) {
                //each mutation has a random position between (1-g)L and (1+g)L
                pop[i].let.push_back(Lv * (1 - gammav + (2 * gammav * rnd.rand())));
                sort(pop[i].let.begin(), pop[i].let.end());
            }

            //// deleterious, same process
            mut = poisdev(Uv * (1 - thv)); // number of new deleterious mutations
            for (j = 0; j < mut; j++) {
                //each mutation has a random position between 0 and 2L
                pop[i].sel.push_back(twoL * rnd.rand());
                sort(pop[i].sel.begin(), pop[i].sel.end());
            }

        } // end of: for (i = 0; i < twoN; i++) { // for each chromosome
		
		
		foutSC.open(nomFichierSC,std::ofstream::app);
		foutSC << SCchm << " ";
		foutSC.close();
		
		
        //// neutrals mutations
		//// same process than before
		
		// linked
		mutne = poisdev(twoN * nbr_locus_neutrev * Unev);
		for(i = 0; i < mutne; i++) {
			chrmut = rnd.randInt(twoN - 1);
			place_chrmut = rnd.randInt(nbr_locus_neutrev - 1);
			neutre = pop[chrmut].neutral_ind[place_chrmut];
			do
            {
                 pop[chrmut].neutral_ind[place_chrmut] = rnd.randInt();
            }
			while (pop[chrmut].neutral_ind[place_chrmut] == neutre);
		}
		// independent
 		mutne = poisdev(twoN * nbr_locus_neutrev * Unev);  
		for(i = 0; i < mutne; i++) { 
			chrmut = rnd.randInt(twoN - 1);
			place_chrmut = rnd.randInt(nbr_locus_neutrev - 1);
			neutre = pop[chrmut].neutral_S[place_chrmut];
			do
            {
                 pop[chrmut].neutral_S[place_chrmut] = rnd.randInt();
            }
			while (pop[chrmut].neutral_S[place_chrmut] == neutre);
		} 

		//// end of mutation phase
		
		//// 
		//// selection, mating, next generation population
		////

		// we need to compute fitness before mating
		// it is now a little bit more complicated because each female part is compatible with some SI pollen (except the two alleles SI it carries) but also with all SC pollen
		// we need to know the exact composition of pollen pool it will depend on fitness, we also need to compute each deme pollen pool
        wbar = 0; // for mean fitness
        wmax = 0; // for maximum fitness
        varw = 0; // for fitness variance

		// for mean fitness by feme
        for (k=0; k < pv ; k++)
            wbar_p[k]=0;

		// for the sum fitness of each SC and SI alleles
        for (un_k = 0; un_k < Pfec.size(); un_k++)
            Pfec[un_k] = 0;

		// for the sum fitness of each SC and SI alleles by deme
        for (int l=0; l < pv; l++) {
            for (un_k = 0; un_k < Pfec.size();un_k++) {
                tab_Pfec[l][un_k]=0;
            }
        }
        // for each individual j: fitness calculation
        for (j = 0; j < Nv; j++) {
            nb = 2 * j; 
            w = fitness(pop[nb], pop[nb + 1], Whet, Whom, Whet_l, Whom_l); // fitness

            Wij[j] = w; // keep safe fitness
            varw += w * w; // for fitness variance

            if (wmax < w) // if it is the highest fitness until now
                wmax = w; // it is now the highest fitness
				
            // sum of fitness for each S allele
            for (k = 0; k < 2; k++) { // for each chromosome of an individual
                Si = pop[nb + k].S; // keep SI allele
                if (Si <0) { // if it is a SC allele
                    Pfec[0] += w; // sum to SC total fitness
                    tab_Pfec[(j-j%taille_patch)/taille_patch][0] += w; // sum to SC deme fitness
                } else { // if it is a SI allele
                    Pfec[Si] += w; // sum to this SI total fitness
                    tab_Pfec[(j-j%taille_patch)/taille_patch][Si] += w; // sum to this SI deme fitness
                }
            }
        } // end of: for (j = 0; j < Nv; j++) { for each individual j: fitness calculation
		
        varw /= Nv; // for fitness variance
        // for fitness mean
        for (un_k = 0; un_k < Pfec.size(); un_k++)
            wbar += Pfec[un_k]; // sum total fitness for each SI allele
		// same by deme
        for (int l=0; l < pv; l++) {
            for (un_k = 0; un_k < Pfec.size();un_k++) {
                wbar_p[l] += tab_Pfec[l][un_k];
            }
        }

		// we have in Pfec the sum of fitness by allele (or the quantity of pollen)
		// we now want the quantity of compatible pollen for each SI allele
		// then we will take the total quantity of pollen and substract the pollen of each SI allele
        Pfec[0] = wbar / 2.0; // for SC, it is compatible with everybody, but we need to divide for some calculation later (Gervais 2014)
                for (un_k = 1; un_k < Pfec.size(); un_k++) // for each SI allele
            Pfec[un_k] = wbar - Pfec[un_k]; // total quantity - quantity of this SI = compatible pollen quantity

		// same process for pollen by deme
        for (int l=0; l < pv; l++) {
            tab_Pfec[l][0]=wbar_p[l]/2;
            for (un_k = 1; un_k < Pfec.size();un_k++) {
                tab_Pfec[l][un_k]=wbar_p[l]-tab_Pfec[l][un_k];
            }
        }
        wbar /= twoN; // we can now get the mean fitness
        varw -= wbar * wbar; // we can finish also for the variance of fitness

        cmptSelf = 0;// for the number of selfing event
        p_disp_bar = 0; // for pollen dispersal
        g_disp_bar = 0; // for seed dispersal
		
        // sampling the next generation:
        for (ind = 0; ind < Nv; ind++) { // for each individual in the next generation population
            ind_chr = 2*ind; // first chromosome of the individual (in even position)
            // sampling the mother:
			// same process than before
            do{ // selection on fitness

                if (rnd.randExc() >= dgv) { // no seed dispersal
                    chrm = (ind_chr-ind_chr%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                    mom = chrm/2; 

                } else { // seed dispersal
                    do{ 
                        chrm = rnd.randInt(twoN-1);
                        mom = chrm/2; 
                    } while (ind_chr-ind_chr%two_taille_patch==chrm-chrm%two_taille_patch); // need an individual from another deme
                }

            } while (rnd.rand() >= Wij[mom] / wmax); // selection on fitness

            // one more seed dispersal (or not)
            if (ind_chr-ind_chr%two_taille_patch!=chrm-chrm%two_taille_patch) {
                g_disp_bar++;
            }

            // recombination for mother gamete:
            if (chrm % 2 == 0)
                rec(temp[ind_chr], pop[chrm], pop[chrm + 1], Lv);
            else
                rec(temp[ind_chr], pop[chrm], pop[chrm - 1], Lv);

            // maternal SI-alleles, as before:
            Si = pop[chrm].S; //
            if(chrm % 2 == 0) {
                Sj = pop[chrm + 1].S;
            } else {
                Sj = pop[chrm - 1].S;
            }
			// if SC alleles, put 0 for all different SC alleles
            if (Sj<0) {
                Sj=0;
            }
            if (Si<0) {
                Si=0;
            }

            //self-incompatibility
            if ((Si == 0) || (Sj == 0)) // if mother have at least on SC allele
				// calculation of the selfing probability for this individual/mother
				// it depend on its pollen production (fitness) and on pollen dispersal rate, pollen pool of the deme and on the metapopulation
                slf = av * Wij[mom] / ( av * Wij[mom] + (1-av) * ( (1-dpv)*((tab_Pfec[(mom-mom%taille_patch)/taille_patch][Si+Sj]-Wij[mom]) / (taille_patch-1) ) + dpv*((Pfec[Si+Sj] - Wij[mom]) / (Nv-1)) ) );
            else // no SC allele for this individual/mother
                slf = 0; // no selfing

            // sampling the father:
            if (rnd.randExc() < slf) { // if it is a selfing case
                cmptSelf++; // one more selfing event

                par1 = 2*mom; // just get chromosomes position

                if ((pop[par1].S < 0) && (pop[par1 + 1].S > 0)) { //if the chromosome in even position have a SC allele but not the one in odd position
                    rec(temp[ind_chr + 1], pop[par1], pop[par1 + 1], Lv); // recombination arround SC allele 
                } else if ((pop[par1].S > 0) && (pop[par1 + 1].S < 0)) { //if chromosome in odd position have a SC allele but not the one in even position
                    rec(temp[ind_chr + 1], pop[par1 + 1], pop[par1], Lv); // recombination arround SC allele 
                } else { // if the individual have two SC alleles
					// randomly choose an allele for recombination
                    if (rnd.rand() < 0.5)
                        rec(temp[ind_chr + 1], pop[par1], pop[par1 + 1], Lv);
                    else 
                        rec(temp[ind_chr + 1], pop[par1 + 1], pop[par1], Lv);
                }

            } else { // it is an outcrossing event, we have the same process than in phase B, but there is some SC pollen in the pool

                do { // until we found a compatible pollen with a sufficient fitness 

                    if (rnd.randExc() >= dpv) {
                        chrd = (chrm-chrm%two_taille_patch)+rnd.randInt(two_taille_patch-1); 
                        dad = chrd/2;
                    } else {
                        do { 
                            chrd = rnd.randInt(twoN-1); 
                            dad = chrd/2; 
                        } while (chrm-chrm%two_taille_patch==chrd-chrd%two_taille_patch);
                    }
				// not the same individual (outcrossing) || fitnes selection || compatibility
                } while ( (mom == dad) || (rnd.rand() >= Wij[dad] / wmax) || (pop[chrd].S == Si && pop[chrd].S * Si != 0) || (pop[chrd].S == Sj && pop[chrd].S * Sj != 0) );
          
			// if there is a pollen dispersal event
            if (chrm-chrm%two_taille_patch!=chrd-chrd%two_taille_patch) {
                p_disp_bar++;
            }
			
            // recombination for father gamete
            if (chrd % 2 == 0)
                rec(temp[ind_chr + 1], pop[chrd], pop[chrd + 1], Lv);
            else
                rec(temp[ind_chr + 1], pop[chrd], pop[chrd - 1], Lv);

            } // end of: if (rnd.randExc() < slf) { // if it is a selfing case } else { it is an outcrossing event

        } // end of: for (ind = 0; ind < Nv; ind++) { // for each individual in the next generation population
		
        // dispersal rate
        g_disp_bar /= Nv;
        p_disp_bar /= Nv;

		// overwrite old metapopulation with the new one
        for (i = 0; i < twoN; i++)
            pop[i] = temp[i];

        //////////////  end of mating, selection, next generation formation
		//////////////***
		//// analysis and results, same process than before with few results on SC
		//////////////

		S_alleles.clear();
		SC_alleles.clear();
		SCch = 0;
        nS_patch_bar=0;
		nS_patch_var=0;
        nSC_patch_bar=0;
		nSC_patch_var=0;
        SC_patch_bar=0;
        SC_patch_var=0;
		for (i = 0; i < twoN; i++) { // for each chromosome
            Si = pop[i].S;
			
            if ( (i%two_taille_patch)==0) { // first chromosome of the deme
                S_alleles_patch.clear(); 
                SC_alleles_patch.clear();
                SCch_patch=0;
            }

            if (Si < 0) { // if it is an SC allele

                SCch_patch++;  // one more SC allele on the deme
				
                for (un_j = 0; un_j < SC_alleles_patch.size(); un_j++)
                    if (Si == SC_alleles_patch[un_j])
                        break; 
                    if (un_j == SC_alleles_patch.size())
                    SC_alleles_patch.push_back(Si);

            } else { // if it is an SI allele
			
                for (un_j = 0; un_j < S_alleles_patch.size(); un_j++)
                    if (Si == S_alleles_patch[un_j]) 
                        break;
                    if (un_j == S_alleles_patch.size())
                    S_alleles_patch.push_back(Si);
            }

            if ( (i%two_taille_patch)==(two_taille_patch-1) ) { // last chromosome of the deme
                nS_patch_bar += S_alleles_patch.size(); // SI allele number
                nS_patch_var += S_alleles_patch.size()*S_alleles_patch.size();

                nSC_patch_bar += SC_alleles_patch.size(); // SC allele number
                nSC_patch_var += SC_alleles_patch.size()*SC_alleles_patch.size();

                SC_patch_bar += SCch_patch/double(two_taille_patch); // SC frequency
                SC_patch_var += (SCch_patch/double(two_taille_patch))*(SCch_patch/double(two_taille_patch));

				// if less than 3 SI alleles and no SC allele
                if (S_alleles_patch.size() < 3 && SCch_patch == 0)
                    signal_d_arret=true;  // stop simulation
            }

            // Stats for the metapopulation, as before
            if (Si < 0) { // SC allele
                SCch++;
                for (un_j = 0; un_j < SC_alleles.size(); un_j++) 
                    if (Si == SC_alleles[un_j]) 
                        break; 
                if (un_j == SC_alleles.size()) 
                    SC_alleles.push_back(Si);
            } else { // SI allele
                for (un_j = 0; un_j < S_alleles.size(); un_j++) 
                    if (Si == S_alleles[un_j]) 
                        break;
                if (un_j == S_alleles.size()) 
                    S_alleles.push_back(Si);
            }
			
        } // end of: for (i = 0; i < twoN; i++) { // for each chromosome
		
		// write specific SC data
		foutSC.open(nomFichierSC,std::ofstream::app);
		foutSC << SCch << " " << SC_alleles.size() << endl;
		foutSC.close();

		// if less than 3 SI alleles and no SC allele
        if (S_alleles.size() < 3 && SCch == 0) 
            signal_d_arret=true; // stop simulation

        //// stats by deme, as before
        nS_patch_bar = nS_patch_bar/double(pv);
        nS_patch_var = (nS_patch_var/double(pv))-nS_patch_bar*nS_patch_bar;
		
        nSC_patch_bar = nSC_patch_bar/double(pv);
        nSC_patch_var = (nSC_patch_var/double(pv))-nSC_patch_bar*nSC_patch_bar;

        SC_patch_bar = SC_patch_bar/double(pv);
        SC_patch_var = (SC_patch_var/double(pv))-(SC_patch_bar*SC_patch_bar);

        // measuring fitness for get inbreeding depression, every "pasv" generations
		// as before
        if (gen % pasv == 0) { 

            le = 0;
            for (i = 0; i < twoN; i++)
                le += pop[i].sel.size();

            le_l = 0;
            for (i = 0; i < twoN; i++)
                le_l += pop[i].let.size();


            if (menagev==1) {
                elim_let = menage_let(pop,twoN);
                elim_del = menage_del(pop,twoN);
            }

			// measuring the average fitness of "sampleSv" selfed offspring,
			// and "sampleSv" offspring produced by outcrossing:
			w_rand = 0; 
			w_out_glo = 0; 
			w_out_loc = 0; 
			w_self = 0; 

			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1);
				chr2 = rnd.randInt(twoN-1); 
                if (chr1 % 2 == 0 ) { 
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv);
				} else { 
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
				}
				if (chr2 % 2 == 0) { 
					rec(off2, pop[chr2], pop[chr2 + 1], Lv);
				} else {
					rec(off2, pop[chr2], pop[chr2 - 1], Lv);
                }
				w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l);
				w_rand += w_m;
				w_rand_var += w_m*w_m;
			}
			w_rand = w_rand/sampleSv;
			w_rand_var = (w_rand_var/sampleSv)-(w_rand*w_rand);

			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1); 
                if (chr1 % 2 == 0 ) {
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv);
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ;
				} else {
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
                    do {
                        chr2 = rnd.randInt(twoN-1);
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) ) ;
				}
				if (chr2 % 2 == 0) {
					rec(off2, pop[chr2], pop[chr2 + 1], Lv);
				} else {
					rec(off2, pop[chr2], pop[chr2 - 1], Lv); 
				}
                w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l); 
				w_out_glo += w_m;
				w_out_glo_var += w_m*w_m;
			}
			w_out_glo = w_out_glo/sampleSv;
			w_out_glo_var = (w_out_glo_var/sampleSv)-(w_out_glo*w_out_glo);

			for (j = 0; j < sampleSv; j++) {
				chr1 = rnd.randInt(twoN-1);
                if (chr1 % 2 == 0 ) {
                    rec(off1, pop[chr1], pop[chr1 + 1], Lv);
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                    } while ( (chr1==chr2) | (chr1+1==chr2) ) ;
				} else {
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
                    do {
                        chr2 = (chr1-chr1%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                    } while ( (chr1==chr2) | (chr1 - 1==chr2) );
				}
				if (chr2 % 2 == 0) 
					rec(off2, pop[chr2], pop[chr2 + 1], Lv);
				else
					rec(off2, pop[chr2], pop[chr2 - 1], Lv);
				w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l);
				w_out_loc += w_m;
				w_out_loc_var += w_m*w_m;
			}
			w_out_loc = w_out_loc/sampleSv;
			w_out_loc_var = (w_out_loc_var/sampleSv)-(w_out_loc*w_out_loc);

			for (j = 0; j < sampleSv; j++) { 
				chr1 = rnd.randInt(twoN-1); 
				if (chr1 % 2 == 0) { 
					nb3 = chr1 + rnd.randInt(1); 
					rec(off1, pop[chr1], pop[chr1 + 1], Lv); 
				} else {
					nb3 = chr1 - rnd.randInt(1);
                    rec(off1, pop[chr1], pop[chr1 - 1], Lv);
				}
				if (nb3 % 2 == 0)
					rec(off2, pop[nb3], pop[nb3 + 1], Lv);
				else
					rec(off2, pop[nb3], pop[nb3 - 1], Lv);
				w_m = fitness(off1, off2, Whet, Whom, Whet_l, Whom_l);
				w_self += w_m;
				w_self_var += w_m*w_m;
			}
			w_self = w_self/sampleSv;
			w_self_var = (w_self_var/sampleSv)-(w_self*w_self);

            /*foutresult <<"generation"
                << " " << "haplotype_number" << " " << "haplotype_diversity" //ligne alleles totales
                << " " << "local_haplotype_number" << " " << "local_haplotype_diversity"  << " " << "local_haplotype_number_variance" << " " << "local_haplotype_diversity_variance" //ligne alleles locale
                << " " << "SC_frequency" << " " << "selfing_rate" << " " << "SC_mean" << " " << "SC_variance" << " " << "SC_haplotype_number" << " " << "SC_haplotype_diversity" << " " << "SC_haplotype_variance"  //ligne alleles SC
                << " " << "fitness_mean" << " " << "fitness_variance" << " " << "mut_per_chr_del" << " " << "mut_per_chr_let" << " " << "fixe_del" << " " << "fix_let" // ligne fitness et mutations
                << " " << "fit_rand" << " " << "fit_out_glo" << " " << "fit_out_glo_SI" << " " << "fit_out_loc" << " " << "fit_out_loc_SI" << " " << "fit_self" // ligne fitness théoriques
                << " " << "var_fit_rand" << " " << "var_fit_out_glo" << " " << "var_fit_out_glo_SI" << " " << "var_fit_out_loc" << " " << "var_fit_out_loc_SI" << " " << "var_fit_self" // ligne var fitness théoriques
                << " " << "rejection_glo" << " " << "rejection_loc" << " " << "rejection_fit_glo" << " " << "rejection_fit_loc" // ligne proba rejet
                << " " << "disp_pol" << " " << "disp_graine" <<  endl; // ligne disp
            foutresult.close();*/

            foutresult.open(nomFichierResult,std::ofstream::app);
			foutresult << NbPrelimAv+NbPrelimBv+gen
                << " " << S_alleles.size() << " " << D //ligne alleles totales
                << " " << nS_patch_bar << " " << D_patch << " " << nS_patch_var << " " << D_patch_var //ligne alleles locale
                << " " << double (SCch) / twoN << " " << double (cmptSelf) / Nv << " " << SC_patch_bar << " " << SC_patch_var << " " << nSC_patch_bar << " " << "NA" << " " << nSC_patch_var   //ligne alleles SC
                << " " << wbar << " " << varw << " " << le / twoN << " " << le_l / twoN << " " << elim_del << " " << elim_let // ligne fitness et mutations
                << " " << w_rand << " "  << w_out_glo << " "  << "NA" << " "  << w_out_loc << " "  << "NA" << " "  << w_self // ligne fitness théoriques
                << " " << w_rand_var << " "  << w_out_glo_var << " "  << w_out_glo_SI_var << " "  << w_out_loc_var << " "  << w_out_loc_SI_var << " "  << w_self_var // ligne var fitness théoriques
                << " " << reject_glo/sampleSv << " " << reject_loc/sampleSv << " " << reject_fit_glo/sampleSv << " " << reject_fit_loc/sampleSv // ligne proba rejet
                << " " << p_disp_bar << " " << g_disp_bar << endl; // ligne disp
            foutresult.close();
			
			// note SC frequency for now if SC have invade.
            if (gen >= NbGenv - 100 * pasv) 
                fSC.push_back(double (SCch) / twoN);

            // neutrals alleles, same than in phase A and B
			foutFST.open(nomFichierFST,std::ofstream::app);
            foutFST << gen;
            for (i = 0; i<nbr_locus_neutrev; i++) {
                identite_lie_glo = 0;
                identite_lie_loc = 0;
                identite_lie_indiv = 0;
                identite_ind_glo = 0;
                identite_ind_loc = 0;
                identite_ind_indiv = 0;
                for (j = 0; j < sampleSv; j++) {
                    chr_foc = rnd.randInt(twoN-1);
                    if (chr_foc % 2 == 0 ) {

                        chr_comp = chr_foc+1;
                        do {
                            chr_glo = rnd.randInt(twoN-1);
                        } while ( (chr_glo==chr_foc) | (chr_glo==chr_comp) ) ;
                        do {
                            chr_loc = (chr_foc-chr_foc%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                        } while ( (chr_loc==chr_foc) | (chr_loc==chr_comp) ) ;

                    } else {

                        chr_comp = chr_foc-1;
                        do { 
                            chr_glo = rnd.randInt(twoN-1);
                        } while ( (chr_glo==chr_foc) | (chr_glo==chr_comp) ) ;
                        do { 
                            chr_loc = (chr_foc-chr_foc%two_taille_patch)+rnd.randInt(two_taille_patch-1);
                        } while ( (chr_loc==chr_foc) | (chr_loc==chr_comp) ) ;

                    }

                    if ( pop[chr_foc].neutral_S[i]==pop[chr_comp].neutral_S[i] ) {
                        identite_lie_indiv++;
                    }
                    if ( pop[chr_foc].neutral_S[i]==pop[chr_glo].neutral_S[i] ) {
                        identite_lie_glo++;
                    }
                    if ( pop[chr_foc].neutral_S[i]==pop[chr_loc].neutral_S[i] ) {
                        identite_lie_loc++;
                    }
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_comp].neutral_ind[i] ) {
                        identite_ind_indiv++;
                    }
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_glo].neutral_ind[i] ) {
                        identite_ind_glo++;
                    }
                    if ( pop[chr_foc].neutral_ind[i]==pop[chr_loc].neutral_ind[i] ) {
                        identite_ind_loc++;
                    }
                } 
				foutFST << " " << identite_lie_indiv << " " << identite_lie_glo << " " << identite_lie_loc << " " << identite_ind_indiv << " " << identite_ind_glo << " " << identite_ind_loc;
            }
            foutFST << endl;
			foutFST.close();

        } // end of: if (gen % pasv == 0) { if it is a generation for write results

		// if there is only SC in the metapopulation
		if (SCch == twoN) 
			cmptgen++; // count it

		// if only SC for sufficient time, we can stop the simulation (no SC->SI so invasion whatever happen
		if (cmptgen >= pasv*11 && cmptgen >= NbGenv/10  ) { /
			fSC.clear();
			fSC.push_back(double (SCch) / twoN); // SC frequency for the "main", if we stop we say it was a full invasion so it is 1
			signal_d_arret=true;
			foutsuivi << "invasion ";
		}
		
		// if we need to stop
		if ( signal_d_arret==true ) { 
			break; 
			foutsuivi << "arret ";
		}

	} // end of:  for ( gen = 0; gen<NbGenv; gen++) { // for each generation until the end

	
	// snapshot at the end of the simulation
    for (i = 0; i < Nv; i++) {
         fitness_temp[i] = fitness(pop[i], pop[i+1], Whet, Whom, Whet_l, Whom_l);
    }
    foutPhoto.open(nomFichierPhoto,std::ofstream::app);
    for (i = 0; i < twoN; i++) {
        foutPhoto << "C" << " " << ((i-(i%two_taille_patch))/two_taille_patch)+1  << " " << pop[i].S << " " << fitness_temp[i/2] << endl;
    }
    foutPhoto.close();

	delete [] pop;
	delete [] temp;
	delete [] Wij;
	fSCbar = 0;
	depbar = 0;

    // mean SC frequency
	for (un_i = 0; un_i < fSC.size(); un_i++)
    fSCbar += fSC[un_i];
    fSCbar /= fSC.size();
    // mean inbreeding depression
    for (un_i = 0; un_i < dep.size(); un_i++)
    depbar += dep[un_i];
    depbar /= dep.size();
    Res.frSC = fSCbar;
    Res.delta = depbar;

	foutsuivi << "fin" << endl;
	foutsuivi.close();

	// return mean SC frequency and mean inbreeding depression to the "main"
    return Res;
}
