// Functions for fitness computation and for recombination:

#include "depression.h"
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

extern MTRand rnd;


//Computes the fitness of a diploid with chromosomes c1 and c2.
//Multiplicative model, wHe = 1-hs, wHo = 1-s:

double fitness(chr &c1, chr &c2, double wHe, double wHo, double wHe_l, double wHo_l) {
// receive: chr 1 and chr 2 from focal individual, fitness cost values 

	int s1 = c1.sel.size() - 1; //number of deleterious mutations on chr1 -1
	int s2 = c2.sel.size() - 1; //number of deleterious mutations on chr2 -1
	int s1_l = c1.let.size() - 1; //number of lethal mutations on chr1 -1
	int s2_l = c2.let.size() - 1; //number of lethal mutations on chr2 -1

	double w = 1.0; //initial fitness = 1

	// it is important to read the entire loop to well understand the role of s1 and s2
	while ((s1 > -1) || (s2 > -1)) { //while deleterious mutations counters are not to -1

		if (s1 == -1) { //if there is no more deleterious on chr1
			w *= pow(wHe, s2 + 1); //fitness decreas du to the remaining heterozygous deleterious mutations on chr2
			break; //going out of the while loop, no more deleterious mutations
		}

		if (s2 == -1) { //if there is no more deleterious on chr2
			w *= pow(wHe, s1 + 1); //fitness decreas du to the remaining heterozygous deleterious mutations on chr1
			break; //going out of the while loop, no more deleterious mutations
		}

		// heterozygous mutation on c1:
		if (c1.sel[s1] > c2.sel[s2]) { //if the position of mut on chr1 is higher than on chr2
			w *= wHe; //fitness decrease
			s1--; //we have "used" one mutation on chr1
			continue; //go back to begining of the while
		}

		// heterozygous mutation on c2:
		if (c1.sel[s1] < c2.sel[s2]) { //if the position of mut on chr2 is higher than on chr1
			w *= wHe; //fitness decrease
			s2--; //we have "used" one mutation on chr2
			continue; //go back to begining of the while
		}

		// homozygous mutation: if any of the previous "if" was true: there is an homozygous mutation
		w *= wHo; //fitness decrease
		s1--; //we have "used" one mutation on chr1
		s2--; //we have "used" one mutation on chr2
	}

	// lethal mutations
	// same process than for deleterious mutations but with wHe_l and wHo_l for fitness decrease
	while ((s1_l > -1) || (s2_l > -1)) { 

		if (s1_l == -1) { 
			w *= pow(wHe_l, s2_l + 1); 
			break; 
		}

		if (s2_l == -1)	{ 
			w *= pow(wHe_l, s1_l + 1); 
			break; 
		}

		// heterozygous mutation on c1:

		if (c1.let[s1_l] > c2.let[s2_l]) {
			w *= wHe_l; 
			s1_l--; 
			continue;
		}

		// heterozygous mutation on c2:

		if (c1.let[s1_l] < c2.let[s2_l]) { 
			w *= wHe_l; 
			s2_l--; 
			continue;
		}

		// homozygous mutation:
		w *= wHo_l;
		s1_l--;
		s2_l--; 
	}

	return w; //on renvoit la valeur de fitness
} //fin de la fonction

// Search for fixed mutation and remove them from the population
// We used the first chromosome of the population as reference, if a mutation is fixed it got it 
int menage_let(chr * population,int taille_pop_chr) {
    int elimination=0; 
    bool nonfixe; // variable that become "true" if the mutation is not fixed
	int compt_chr; // counter for each chromosome
    bool trouve; // become true if we found the focal mutation on the focal chromosome
	unsigned int place; // position of the mutation on chromosomes
    int note[taille_pop_chr]; // to keep position of each mutation to be eliminated
    // for each mutation of the first chromosome of the first individual of the population
    for (unsigned int i=0; i < population[0].let.size(); i++) {
        // keep the position of the focal mutation 
        double mut_focale = population[0].let[i];
        note[0]=i;
        nonfixe=false; // we don t know if the mutation is not fixed so it is false
        compt_chr = 0; // we will start to investigate the chromosome at position 0+1
		// while we dont have a proof that the mutation is not fixed and that it remains some individuals to investigate
        while ( (nonfixe==false) & (compt_chr < taille_pop_chr-1) ) {
            compt_chr++; // next chromosome 
            trouve = false; // initialisation (we still did not find the mutation for this chromosome)
            place = 0; // initialisation of mutation position
            // while we dont have find the mutation and that we haven't gone beyond the mutation position
            while ( (trouve==false) & (place < population[compt_chr].let.size()) ) {
                //  if it is the same mutation
                if ( population[compt_chr].let[place] == mut_focale ) {
                    // we found the mutation
                    trouve=true;
					// keep the identity
                    note[compt_chr]=place;
                }
                // on incrémente le compteur qui parcoure le chr
                place++;
            }
            // if we didn't find the focal mutation
            if (trouve==false) {
                nonfixe=true; // non fixed mutation
            }
        } 
        // if fixed mutation
        if (nonfixe==false) {
            // for each chromosome of the population
            for (int j = 0; j < taille_pop_chr; j++) {
                // erase the mutation
                population[j].let.erase(population[j].let.begin()+note[j]);
            } 
            elimination++; // we had erased one mutation
        } 
    } 
    return(elimination);
} // function end

// exactly the same process than for menage_let just above
int menage_del(chr * population,int taille_pop_chr) {
    int elimination=0;
    bool nonfixe; 
    int compt_chr;
    bool trouve;
    unsigned int place;
    int note[taille_pop_chr];
    for (unsigned int i=0; i < population[0].sel.size(); i++) {
        double mut_focale = population[0].sel[i];
        note[0]=i;
        nonfixe=false;
        compt_chr = 0;
        while ( (nonfixe==false) & (compt_chr < taille_pop_chr-1) ) { 
            compt_chr++;
            trouve = false;
            place = 0; 
            while ( (trouve==false) & (place < population[compt_chr].sel.size()) ) {
                if ( population[compt_chr].sel[place] == mut_focale ) {
                    trouve=true;
                    note[compt_chr] = place;
                }
                place++;
            }
            if (trouve==false) {
                nonfixe=true;
            }
        }
        if (nonfixe==false) {
            for (int j = 0; j < taille_pop_chr; j++) {
                population[j].sel.erase(population[j].sel.begin()+note[j]);
            }
        elimination++;
        }
    }
    return(elimination);
}

// Constructs chromosome "res" by recombining chromosomes c1 and c2
// "Sz" is the map length L (average number of cross-overs per meiosis).
// The S-locus is inherited from chromosome c1.
// The number and position of cross-overs are sampled randomly.
void rec(chr &res, chr &c1, chr &c2, double Sz) { 
// receive temporary place , chr1, chr2 and the number of crossover

	int sz1 = c1.sel.size(); //number of del mutations chr1
	int sz2 = c2.sel.size(); //number of del mutations chr2
	int s1 = sz1 - 1; 
	int s2 = sz2 - 1; 

	int sz1_l = c1.let.size(); //number of let mutations chr1
	int sz2_l = c2.let.size(); //number of let mutations chr2
	int s1_l = sz1_l - 1;
	int s2_l = sz2_l - 1;

	res.S = c1.S; // haplotype S inheritance (always chr1)

    // sub fonction for neutral alleles 
	// neutral_S are totaly linked to S locus
    res.neutral_S.clear(); // clear
    res.neutral_S = c1.neutral_S; // inherite neutral from chr1 (chr1 inherite the haplotype S)

	// neutral_ind are totaly independend of S locus
    res.neutral_ind.clear();
    for (unsigned int l=0; l < c1.neutral_ind.size();l++) {
        if (rnd.randInt(1)==0) {
            res.neutral_ind.push_back(c1.neutral_ind[l]);
        } else {
            res.neutral_ind.push_back(c2.neutral_ind[l]);
        }
    }

    // clear vectors
    res.sel.clear();
    res.let.clear();

    int j; //
    chr Co; //

    double twoS = 2 * Sz; // twice the average number of crossing over/chromosome size

    // number of cross-overs:
    int nbCo = int(poisdev(Sz));

    // positions of cross-overs (between 0 and 2L) are put in the vector Co.sel
    for (j = 0; j < nbCo; j++) //for each crossing over
    Co.sel.push_back(twoS * rnd.rand()); // keep crossing over position

    sort(Co.sel.begin(), Co.sel.end()); // order crossing over positions

    // Counting the number of cross-overs on the right of the S-locus:
    int cmpt = 0; // counter set up
    j = Co.sel.size() - 1; //number/counter of crossing over - 1

    while ((j > -1) && (Co.sel[j] > Sz)) { //while counter of crossing over > -1 and position > Sz (midle of the chromosome)
        cmpt++; //a grossing over at the "right" of the chromosome
        j--;// a crossing over was done
    }

    // if the number of cross-overs on the right is even:
    if (cmpt % 2 == 0) { 

        // for each cross-over (starting on extreme right)
        for (j = 1; j <= cmpt; j++)	{ //for each corssing over
            if (j % 2 == 1)	{ //if its odd number

                // all mutations on the right of cross-over on chromosome 1 are incorporated

                // for deleterious mutations
                while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j])) { // while there is deleterious mutation & we are still at the right of the chromosome
                    res.sel.push_back(c1.sel[s1]); //one more deleterious for the new chromosome
                    s1--;
                }
                while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j])) // while there is deleterious mutation & we are still at the right of the chromosome
                    s2--; //we need to discard these mutations too

                // for lethal mutation, same process
                while ((s1_l > -1) && (c1.let[s1_l] > Co.sel[nbCo - j])) {
                    res.let.push_back(c1.let[s1_l]);
                    s1_l--;
                }
                while ((s2_l > -1) && (c2.let[s2_l] > Co.sel[nbCo - j]))
                    s2_l--;

            } else { //j is even number
                // all mutations on the right of cross-over on chromosome 2 are incorporated

                // for deleterious mutations, same process than before
                while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j])) {
                    res.sel.push_back(c2.sel[s2]);
                    s2--;
                }
                while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j]))
                    s1--;

                // for lethal mutation, same process
                while ((s2_l > -1) && (c2.let[s2_l] > Co.sel[nbCo - j])) {
                    res.let.push_back(c2.let[s2_l]);
                    s2_l--;
                }
                while ((s1_l > -1) && (c1.let[s1_l] > Co.sel[nbCo - j]))
                    s1_l--;

            }

        }

    } else {  // if the number of cross-overs on the right is odd:
	// same general process than before

        // for each cross-over (starting on extreme right):
        for (j = 1; j <= cmpt; j++) {
            if (j % 2 == 1) {
                // all mutations on the right of cross-over on chromosome 2 are incorporated

                // deleterious mutations
                while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j])) {
                    res.sel.push_back(c2.sel[s2]);
                    s2--;
                }
                while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j]))
                    s1--;

                // lethals mutations
                while ((s2_l > -1) && (c2.let[s2_l] > Co.sel[nbCo - j])) {
                    res.let.push_back(c2.let[s2_l]);
                    s2_l--;
                }
                while ((s1_l > -1) && (c1.let[s1_l] > Co.sel[nbCo - j]))
                    s1_l--;


            } else {
                // all mutations on the right of cross-over on chromosome 1 are incorporated

                // deleterious mutations
                while ((s1 > -1) && (c1.sel[s1] > Co.sel[nbCo - j])) {
                    res.sel.push_back(c1.sel[s1]);
                    s1--;
                }
                while ((s2 > -1) && (c2.sel[s2] > Co.sel[nbCo - j]))
                    s2--;

                // lethals mutations
                while ((s1_l > -1) && (c1.let[s1_l] > Co.sel[nbCo - j])) {
                    res.let.push_back(c1.let[s1_l]);
                    s1_l--;
                }
                while ((s2_l > -1) && (c2.let[s2_l] > Co.sel[nbCo - j]))
                    s2_l--;

            } // end of: if (j % 2 == 1)

        } // end of: for (j = 1; j <= cmpt; j++) {

    } // end of: if (cmpt % 2 == 0) { 


    // mutations between the S-locus and the nearest cross-over on the right (on chromosome 1):
    while ((s1 > -1) && (c1.sel[s1] > Sz)) {
        res.sel.push_back(c1.sel[s1]);
        s1--;
    }
    while ((s2 > -1) && (c2.sel[s2] > Sz))
        s2--;

    while ((s1_l > -1) && (c1.let[s1_l] > Sz))	{
        res.let.push_back(c1.let[s1_l]);
        s1_l--;
    }
    while ((s2_l > -1) && (c2.let[s2_l] > Sz))
        s2_l--;

    // number of cross-overs on the left of the S-locus:
    int frst = nbCo - cmpt;

	// same process than before for the left part
    // for each cross-over (starting with the nearest to the S-locus on the left):
    for (j = 1; j <= frst; j++) {
        if (j % 2 == 1) {
            // all mutations on the right of cross-over on chromosome 1 are incorporated
            while ((s1 > -1) && (c1.sel[s1] > Co.sel[frst - j])) {
                res.sel.push_back(c1.sel[s1]);
                s1--;
            }
            while ((s2 > -1) && (c2.sel[s2] > Co.sel[frst - j]))
                s2--;

            while ((s1_l > -1) && (c1.let[s1_l] > Co.sel[frst - j])) {
                res.let.push_back(c1.let[s1_l]);
                s1_l--;
            }
            while ((s2_l > -1) && (c2.let[s2_l] > Co.sel[frst - j]))
                s2_l--;

        } else {
            // all mutations on the right of cross-over on chromosome 2 are incorporated

            while ((s2 > -1) && (c2.sel[s2] > Co.sel[frst - j])) {
                res.sel.push_back(c2.sel[s2]);
                s2--;
            }
            while ((s1 > -1) && (c1.sel[s1] > Co.sel[frst - j]))
                s1--;

            while ((s2_l > -1) && (c2.let[s2_l] > Co.sel[frst - j])) {
                res.let.push_back(c2.let[s2_l]);
                s2_l--;
            }
            while ((s1_l > -1) && (c1.let[s1_l] > Co.sel[frst - j]))
                s1_l--;

        }
    }

    // mutations on the left of the left-most cross-over:
    if (frst % 2 == 0) {

        while (s1 > -1) {
            res.sel.push_back(c1.sel[s1]);
            s1--;
        }

        while (s1_l > -1) {
            res.let.push_back(c1.let[s1_l]);
            s1_l--;
        }

    } else {

        while (s2 > -1) {
            res.sel.push_back(c2.sel[s2]);
            s2--;
        }

        while (s2_l > -1) {
            res.let.push_back(c2.let[s2_l]);
            s2_l--;
        }

    }

	// sorts mutations on offspring chromosome:
	sort(res.sel.begin(), res.sel.end());
	sort(res.let.begin(), res.let.end());

}


///

