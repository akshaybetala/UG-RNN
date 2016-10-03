#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "General.h"
#include <stack>
#include <queue>
#include <string>
#include "Molecule.h"

		Molecule::Molecule(istream &is) {

			is >> foo_string;
			//@agrume			
			//cout << "hello agrume " << "foo_string is:" << foo_string << endl;
			//cout << foo_string << "\n";
			/////////
			is >> foo_string; /* Not always there */
			//@agrume						
			//cout << "hello agrume " << "foo_string is:" << foo_string << endl;
										
			is >> atoms;
			//@agrume 
			//cout << "hello agrume atoms=" <<atoms << endl;
			is >> bonds;
			//cout << "hello agrume bonds=" << bonds << endl;
			is >> foo_string;
			//cout << "hello agrume foo_string=" << foo_string << endl;
			while (foo_string.compare("V2000") != 0) {
				is >> foo_string;
			}

			if (atoms==0 && bonds == 0) {
				cout << "No Atoms! No Bonds! what daaa ffaaaa\n";
				exit(0);
			}

			//@agrume: set the properties for each atom in the molecule
			x = new Float[atoms];
			memset(x, 0, sizeof(Float)*atoms);
			y = new Float[atoms];
			memset(y, 0, sizeof(Float)*atoms);
			z = new Float[atoms];
			memset(z, 0, sizeof(Float)*atoms);
			u = new int[atoms];
			memset(u, 0, sizeof(int)*atoms);
			iso = new int[atoms];
			memset(iso, 0, sizeof(int)*atoms);			
			c = new int[atoms];
			memset(c, 0, sizeof(int)*atoms);
			r = new int[atoms];
			memset(r, 0, sizeof(int)*atoms);
			stereo = new int[atoms];
			memset(stereo, 0, sizeof(int)*atoms);
			valancy = new int[atoms];
			memset(valancy, 0, sizeof(int)*atoms);
			hs = new int[atoms];
			memset(hs, 0, sizeof(int)*atoms);
			aroma = new int[atoms];
			memset(aroma, 0, sizeof(int)*atoms);
			anti = new int[atoms];
			memset(anti, 0, sizeof(int)*atoms);

			atomchar = new std::string[atoms];
			atomletter = new std::string[atoms];
			/* atoms block */
			std::string ch;
			std::string letter;
			for (int a = 0; a < atoms; a++) {
				is >> x[a] >> y[a] >> z[a];
				//@agrume
				//cout << x[a] << " " << y[a] << "\n";
				is >> ch;
				//@agrume				
				//cout << ch << endl;				
				//@agrume
				u[a] = translateU(ch);
				
				/* atom letters */
				letter = "";
				for(int l = 0; (l < ch.length() && ch[l]!='.'); l++)
					letter += ch[l];
				//@agrume				
				//cout << letter << "\n";
				atomchar[a] = letter;
				atomletter[a] = "";
				for(int l = 0; (l < ch.length()); l++) 
					if (ch[l]!='.')
						atomletter[a] += ch[l];
						
						
				valancy[a] = getValance(ch);
				is >> iso[a];
				is >> c[a];
				//@agrume what does it mean?
				if (c[a] == 1) 
					c[a] = 3;
				if (c[a] == 2) 
					c[a] = 2;
				if (c[a] == 3) 
					c[a] = 1;
				if (c[a] == 4) 
					c[a] = 0;
				if (c[a] == 5) 
					c[a] = -1;
				if (c[a] == 6) 
					c[a] = -2;
				if (c[a] == 7) 
					c[a] = -3;

				is >> stereo[a];
				if (stereo[a] > 0) 
					cout << "STEREO HERE *************\n";

				for (i = 0; i < 2; i++)
					is >> foo_int;
			}			
			
			/* bonds block */
			b_type = new int*[atoms];
			b_stereo = new int*[atoms];
			for (i = 0; i < atoms; i++) {
				b_type[i] = new int[atoms];
				memset(b_type[i], 0, sizeof(int)*(atoms));
				//@agrume why atoms+1?					
				b_stereo[i] = new int[(atoms+1)];
				memset(b_stereo[i], 0, sizeof(int)*(atoms));
			}
			
			int j,k;				
			for (i = 0; i < bonds; i++) {
				is >> j;
				//@agrume j is the number of atom
				//cout << "j=" << j <<endl;				
				is >> k;
				//@agrume k is the number of atom to which atom j is bounded
				//cout << "k=" << k <<endl;
				is >> b_type[j-1][k-1];
				//@agrume bond type: 1=single_bond , 2=covalent_bond				
				if(b_type[j-1][k-1]==0) b_type[j-1][k-1]=1;
				b_type[k-1][j-1] = b_type[j-1][k-1];
				//cout << j << " "<< " " << k << " " << b_type[j-1][k-1] << endl;
				//@agrume information about stereo-isomery 				
				is >> b_stereo[j-1][k-1];
				b_stereo[k-1][j-1] = b_stereo[j-1][k-1];
				for (int m = 0; m < 2; m++) 
					is >> foo_int;
			}			
			/* properties */
			is >> foo_string;
			//@agrume
			//cout << "foo_string=" << foo_string << endl;
			int num_chg;
			int num_rads;
			int num_iso;
			while (foo_string.compare("END") != 0) {
				if (foo_string.compare("CHG") == 0) {
					is >> num_chg;
					for(i = 0; i < num_chg; i++) {
						is >> foo_int;
						if (c[foo_int-1] > 0) 
							is >> c[foo_int-1];
					}
					is >> foo_string;
				}
				else if (foo_string.compare("RAD") == 0) {
					cout << "RADICALS ********** \n";
					is >> num_rads;
					for(i = 0; i < num_rads; i++) {
						is >> foo_int;					
						is >> r[foo_int-1];
					}
					is >> foo_string;
					is >> foo_string;
				}
				else if (foo_string.compare("ISO") == 0) {
					cout << "ISOTOPES ****** \n";
					is >> num_iso;
					for(i = 0; i < num_iso; i++) {
						is >> foo_int;
						if (iso[foo_int] > 0) 
							is >> iso[foo_int-1];
					}
					is >> foo_string;					
				}
				else if (foo_string.compare("M") == 0 ){
					is >> foo_string;
				}
				else if (foo_string.compare("A") == 0) {
					//is >> foo_int;
					//is >> foo_string;
					//is >> foo_string;
				}
				else {
					cout << "Unknown property " << foo_string << "\n";
					exit(0);
				}
			}
			
			is >> foo_string;
			//is >> foo_string; is >> foo_string; is >> foo_string;
			//@agrume what is target?
			is >> target;
			//cout << "target" << target << "\n";
			is >> logP;
			//cout << "logP=" << logP << endl;
			
			
			//cout << "Creating nodes \n" << flush;
			all_atoms = new Node*[atoms];
			for (int a = 0; a < atoms; a++) {
				all_atoms[a] = new Node(a, u[a], iso[a], c[a], r[a], stereo[a], valancy[a], x[a], y[a], z[a], aroma[a]);
				//all_atoms[a]->setHydrogen(hs[a]);
			}
			
			//cout << "Creating Graph \n" << flush;
			g = new Graph(all_atoms, b_type, b_stereo, atoms, bonds);			
			//g->orderGraph();	
			//g->loadNeighbours();
			//g->lightenGraph();	
			/* this locates the neighbourhood a particular atom belongs to */
			atomMap = new int[atoms];
			//@agrume what does it do?			
			for (int a = 0; a < atoms; a++)
				atomMap[a] = -1;
			//cout << "Delations \n" << flush;
			//@agrume is it a joke? 
			//@agrume I think they are not deleted in graph...
			/* delete all temp arrays */
			//for (int a = 0; a < atoms; a++) {
			//	//delete all_atoms[a];	/* they get deleted in Graph im not too sure about this */
			//	//delete[] b_type[a];
			//	//delete[] b_stereo[a];
			//}
			//delete all_atoms;
			//delete[] b_type;
			//delete[] b_stereo;
			delete[] aroma;
			delete[] u;
			delete[] iso;
			delete[] c;
			delete[] r;
			delete[] stereo;
			delete[] x;
			delete[] y;
			delete[] z;
			delete[] hs;
		//@agrume here finish the constructor for molecule			
		}
Molecule::~Molecule(){
	//@debug
        //cout << "calling Molecule destructor" << endl;
        ////////
	for(int a=0;a<atoms;a++){
		delete[] b_type[a];
		delete[] b_stereo[a];
		delete all_atoms[a];
	}
	delete[] all_atoms;
	delete[] b_type;
	delete[] b_stereo;
	delete g;		
}

void Molecule::lightenMolecule(){
	delete[] atomchar;
	delete[] atomletter;
	delete[] atomMap;
	delete[] valancy;
	//for(int a=0;a<atoms;a++){
	//	delete[] b_type[a];
	//	delete[] b_stereo[a];
	//	delete all_atoms[a];
	//}
	//delete[] all_atoms;
	//delete[] b_type;
	//delete[] b_stereo;
}

void Molecule::loadweights() {
	ifstream is("avgweights");

	weights = new Float[atoms];

	for (int i = 0; i < atoms; i++)
		is >> weights[i];

	Float max = 0;
	for (int i = 0; i < atoms; i++) {
		if (weights[i] > max) {
			max = weights[i];
		}
	}

	for (int i = 0; i < atoms; i++) {
//		weights[i] = weights[i]/(double)max;
		weights[i] = 1.0;
	}
}

void Molecule::getMol2Info(istream &is) {
			mol2Charges = new Float[atoms];
			memset(mol2Charges, 0, sizeof(Float)*atoms);

			is >> foo_int;
			is >> foo_int;
			is >> foo_int;
			is >> foo_string;
			is >> foo_string;
			is >> foo_string;
			is >> foo_string;
			is >> foo_int;		/* energy value */
			is >> foo_string;
			std::string ch;
			formula = "";
			std::string letter;
			for (int a = 0; a < atoms; a++) {
				is >> foo_int;
				is >> foo_string;
				is >> x[a] >> y[a] >> z[a];
				is >> ch;
				//atomchar[a] = ch;
//				cout << ch << " ";
				formula += ch + " ";
				u[a] = translateMol2U(ch);
				letter = "";
				for(int l = 0; (l < ch.length() && ch[l]!='.'); l++)
					letter += ch[l];
				//cout << letter << "\n";
				atomchar[a] = letter;
				atomletter[a] = "";
				for(int l = 0; (l < ch.length()); l++) 
					if (ch[l]!='.')
						atomletter[a] += ch[l];
				

				atomletter[a]=atomchar[a];
				valancy[a] = getValance(letter);
				is >> foo_int;
				is >> foo_string;
				is >> mol2Charges[a];
			}
//			cout << "\n";
			is >> foo_string;
			int j;
			int k;
			for (i = 0; i < bonds; i++) {
				is >> j;
				is >> j;
				is >> k;
				is >> ch;
				if (ch.compare("1")==0)
					b_type[j-1][k-1] = 1;
				else if (ch.compare("2")==0)
					b_type[j-1][k-1] = 2;
				else if (ch.compare("3")==0)
					b_type[j-1][k-1] = 3;
				else if (ch.compare("ar")==0)
					b_type[j-1][k-1] = 4;
				else if (ch.compare("am")==0)
					b_type[j-1][k-1] = 5;
				else 
					b_type[j-1][k-1] = 6;

				b_type[k-1][j-1] = b_type[j-1][k-1];
			}
}

int Molecule::translateMol2U(std::string s) {
	if (s.compare("C.3") == 0){ 
	return 0;
	} 
	else if (s.compare("C.2") == 0){
	return 1; 
	} 
	else if (s.compare("C.1") == 0){	
	return 2;
	} 
	else if (s.compare("C.ar") == 0){	
	return 3;
	} 
	else if (s.compare("C.cat") == 0){	
	return 4;
	} 
	else if (s.compare("N.3") == 0){	
	return 5;
	} 
	else if (s.compare("N.2") == 0){	
	return 6;
	} 
	else if (s.compare("N.1") == 0){	
	return 7;
	} 
	else if (s.compare("N.ar") == 0){	
	return 8;
	} 
	else if (s.compare("N.am") == 0){
	return 9;
	} 
	else if (s.compare("N.pl3") == 0){	
	return 10;
	} 
	else if (s.compare("N.4") == 0){	
	return 11;
	} 
	else if (s.compare("O.3") == 0){	
	return 12;
	} 
	else if (s.compare("O.2") == 0){	
	return 13;
	} 
	else if (s.compare("O.co2") == 0){	
	return 14;
	} 
	else if (s.compare("O.spc") == 0){	
	return 15;
	} 
	else if (s.compare("O.tp3") == 0){	
	return 16;
	} 
	else if (s.compare("S.3") == 0){	
	return 17;
	} 
	else if (s.compare("S.2") == 0){	
	return 18;
	} 
	else if (s.compare("S.O") == 0){ 	
	return 19;
	} 
	else if (s.compare("S.o") == 0){ 	
	return 20;
	} 
	else if (s.compare("S.O2") == 0){	
	return 21;
	} 
	else if (s.compare("S.o2") == 0){	
	return 22;
	} 
	else if (s.compare("P") == 0){	
	return 23;
	} 
	else if (s.compare("P.3") == 0){	
	return 24;
	} 
	else if (s.compare("F") == 0){	
	return 25;
	} 
	else if (s.compare("Cl") == 0){	
	return 26;
	} 
	else if (s.compare("Br") == 0){	
	return 27;
	} 
	else if (s.compare("I") == 0){	
	return 28;
	} 
	else if (s.compare("In") == 0){	
	return 28;
	} 
	else if (s.compare("H") == 0){	
	return 56;
	} 
	else if (s.compare("H.spc") == 0){	
	return 56;
	} 
	else if (s.compare("H.t3p") == 0){	
	return 56;
	} 
	else if (s.compare("Li") == 0){
	return 29;
	} 
	else if (s.compare("Na") == 0){	
	return 30;
	} 
	else if (s.compare("Mg") == 0){	
	return 31;
	} 
	else if (s.compare("Al") == 0){	
	return 32;
	} 
	else if (s.compare("Si") == 0){	
	return 33;
	} 
	else if (s.compare("K") == 0){	
	return 34;
	} 
	else if (s.compare("Ca") == 0){	
	return 35;
	} 
	else if (s.compare("Cr.th") == 0){	
	return 36;
	} 
	else if (s.compare("Cr.oh") == 0){	
	return 36;
	} 
	else if (s.compare("Mn") == 0){	
	return 37;
	} 
	else if (s.compare("Fe") == 0){
	return 38;
	} 
	else if (s.compare("Co.oh") == 0){	
	return 39;
	} 
	else if (s.compare("Cu") == 0){	
	return 40;
	} 
	else if (s.compare("Zn") == 0){	
	return 41;
	} 
	else if (s.compare("Se") == 0){	
	return 42;
	} 
	else if (s.compare("Mo") == 0){	
	return 43;
	} 
	else if (s.compare("Sn") == 0){	
	return 44;
	} 
	else if (s.compare("Ba") == 0){	
	return 45;
	} 
	else if (s.compare("Sb") == 0){	
	return 46;
	} 
	else if (s.compare("Be") == 0){	
	return 47;
	} 
	else if (s.compare("B") == 0){
	return 48;
	} 
	else if (s.compare("Cd") == 0){	
	return 49;
	} 
	else if (s.compare("As") == 0){	
	return 50;
	} 
	else if (s.compare("Te") == 0){	
	return 51;
	} 
	else if (s.compare("Hg") == 0){	
	return 52;
	} 
	else if (s.compare("Ni") == 0){	
	return 53;
	} 
	else if (s.compare("Ti") == 0){	
	return 54;
	} 
	else if (s.compare("Pb") == 0){	
	return 55;
	} 
	else if (s.compare("LP") == 0){	
	return -1;
	} 
	else if (s.compare("Du") == 0){	
	return -1;
	} 
	else if (s.compare("Du.C") == 0){	
	return -1;
	} 
	else if (s.compare("Any") == 0) {
	return -1;
	} 
	else if (s.compare("Hal") == 0){	
	return -1;
	} 
	else if (s.compare("Het") == 0){	
	return -1;
	} 
	else if (s.compare("Hev") == 0){ 
	return -1;
	} 
	else { 				
		cout << "Unknown symbol " << s << "\n";
		exit(0);
	}
}

		int Molecule::translateU(std::string s) {
			if (s.compare("C") == 0){return 0;} 
			else if (s.compare("H") == 0){return 1;}
			else if (s.compare("N") == 0){return 2;} 
			else if (s.compare("O") == 0){return 3;} 
			else if (s.compare("S") == 0){return 4;} 
			else if (s.compare("F") == 0){return 5;} 
			else if (s.compare("Cl") == 0){return 6;} 
			else if (s.compare("Br") == 0){return 7;}
			else if (s.compare("I") == 0){return 8;} 
			else if (s.compare("In") == 0){return 9;} 			
			else if (s.compare("K") == 0){return 10;}
			else if (s.compare("Na") == 0){return 11;}	//1st salt
			else if (s.compare("Ba") == 0){return 12;}	//2nd salt		
			else if (s.compare("Sb") == 0){return 13;}	//purple
			else if (s.compare("P") == 0){return 14;}	//green
			else if (s.compare("Be") == 0){return 15;}	
			else if (s.compare("Sn") == 0){return 16;}
			else if (s.compare("Cu") == 0){return 17;}	//T metals
			else if (s.compare("B") == 0){return 18;}
			else if (s.compare("Cd") == 0){return 19;}
			else if (s.compare("Ca") == 0){return 20;}
			else if (s.compare("As") == 0){return 21;}
			else if (s.compare("Co") == 0){return 22;}
			else if (s.compare("Cr") == 0){return 23;}
			else if (s.compare("Te") == 0){return 24;}
			else if (s.compare("Fe") == 0){return 25;}
			else if (s.compare("Pb") == 0){return 26;}
			else if (s.compare("Mn") == 0){return 27;}
			else if (s.compare("Hg") == 0){return 28;}
			else if (s.compare("Mo") == 0){return 29;}
			else if (s.compare("Ni") == 0){return 30;}
			else if (s.compare("Se") == 0){return 31;}
			else if (s.compare("Ti") == 0){return 32;}
			else if (s.compare("Zn") == 0){return 33;}
			else if (s.compare("Si") == 0){return 34;}
			else if (s.compare("?") == 0){return 35;}
			else if (s.compare("Mg") == 0){return 36;}
			else if (s.compare("V") == 0){return 37;} 
			else if (s.compare("Li") == 0){return 38;} 
			else if (s.compare("Al") == 0){return 39;}
			else if (s.compare("Zr") == 0){return 40;}
			else if (s.compare("Bi") == 0){return 41;}
			else if (s.compare("Pd") == 0){return 42;}
			else if (s.compare("Pt") == 0){return 43;}
			else if (s.compare("Ru") == 0){return 44;}
			else if (s.compare("Rh") == 0){return 45;}			
			else if (s.compare("Ga") == 0){return 46;}						
			else if (s.compare("Ge") == 0){return 47;}								
			else if (s.compare("Ag") == 0){return 48;}
			else if (s.compare("Tb") == 0){return 49;}
			else if (s.compare("Ir") == 0){return 50;}
			else if (s.compare("W") == 0){return 51;}
			else if (s.compare("Cs") == 0){return 52;}
			else if (s.compare("Re") == 0){return 53;}
			else if (s.compare("Pr") == 0){return 54;}
			else if (s.compare("Nd") == 0){return 55;}
			else if (s.compare("Gd") == 0){return 56;}
			else if (s.compare("Yb") == 0){return 57;}
			else if (s.compare("Er") == 0){return 58;}
			else if (s.compare("U") == 0){return 59;}
			else if (s.compare("Tl") == 0){return 60;}
			else if (s.compare("Au") == 0){return 61;}
			else if (s.compare("Ac") == 0){return 62;}
			else if (s.compare("Ho") == 0){return 63;}
			else if (s.compare("Os") == 0){return 64;}
			else if (s.compare("Sm") == 0){return 65;}
			else if (s.compare("Nb") == 0){return 66;}
			else if (s.compare("R3") == 0){return 67;}
			else if (s.compare("R4") == 0){return 68;}
			else if (s.compare("R5") == 0){return 69;}
			else if (s.compare("R6") == 0){return 70;}
			else if (s.compare("R7") == 0){return 71;}
			else if (s.compare("R8") == 0){return 72;}
			else if (s.compare("R14") == 0){return 73;}
			else if (s.compare("R24") == 0){return 74;}
			else { 					//return 73;
				cout << "***************** WARNING : Unknown symbol " << s << "\n******************\n";
				exit(0);
			}
		}

		/* Determine the cycles from the connection table. Can it be done?*/
		void Molecule::getCycles(int* u, int* k) {
			/* sort bonds here */
			/* j atom must be less than k atom */
			int** atomcon = new int*[bonds]; 
			for (int i = 0; i < bonds; i++) {
				atomcon[i] = new int[10];	/* cant have more than 10 valency */
				for (int j = 0; j < 10; j++) 
					atomcon[i][j] = -1;
			}
			int* val = new int[atoms];
			memset(val, 0, sizeof(int)*atoms);
			for (int i = 0; i < bonds; i++) {
				if (u[i] <= k[i]) {
					atomcon[u[i]-1][val[u[i]]] = k[i]-1;
					val[u[i]]++;
				}
				else {
					atomcon[k[i]-1][val[k[i]]] = u[i]-1;
					val[k[i]]++;				
				}
			}
			/* sort atoms adjacent to current atom */
			int cur, min, min_index, tmp;
			for (int i = 0; i < bonds; i++) {
				for (int n = 0; n < 10; n++) {
					min = 1000;
					min_index=-1;
					for (int j = n; j < 10; j++) {
						if (atomcon[i][j]>=0)
							cur = atomcon[i][j];
						if (cur < min) {
							min = cur;
							min_index = j;
						}
					}		
					tmp = atomcon[i][n];
					atomcon[i][n] = min;
					atomcon[i][min_index] = tmp;
				}
			}
			cout << "The atoms and there connections \n";
			for (int i = 0; i < atoms; i++) {
				for (int j = 0; j < 10; j++) {
					if (atomcon[i][j]>=0)
						cout << i+1 << " " << atomcon[i][j]+1 << "\n";
				}
			}

			int n = bonds;
			n--;
			int cl = 0;
			int looplength = 0;
			int moveto = -1;
			while (n >=0) {
				moveto = u[n];
				while (moveto != k[n] && n>=0) {			
					n--;
				}	
				if (n >= 0) { 
					looplength++;
//					cycles[sc][whichcycle][cl] = u[n];
					cl++;
				}
			}
		}

		int Molecule::getValance(std::string s) {
			if (s.compare("C") == 0){return 4;} 
			else if (s.compare("H") == 0){return 1;}
			else if (s.compare("N") == 0){return 3;} 
			else if (s.compare("O") == 0){return 2;} 
			else if (s.compare("S") == 0){return 2;} 
			else if (s.compare("F") == 0){return 1;} 
			else if (s.compare("Cl") == 0){return 1;} 
			else if (s.compare("Br") == 0){return 1;}
			else if (s.compare("I") == 0){return 1;}
			else if (s.compare("In") == 0){return 1;}  
			else if (s.compare("V") == 0){return 1;}  
			else if (s.compare("K") == 0){return 1;}
			else if (s.compare("Na") == 0){return 1;}
			else if (s.compare("Ba") == 0){return 2;}
			else if (s.compare("Sb") == 0){return 3;}
			else if (s.compare("P") == 0){return 3;}
			else if (s.compare("Be") == 0){return 2;}
			else if (s.compare("Sn") == 0){return 4;}
			else if (s.compare("Cu") == 0){return 3;}
			else if (s.compare("B") == 0){return 3;}
			else if (s.compare("Cd") == 0){return 2;}
			else if (s.compare("Ca") == 0){return 2;}
			else if (s.compare("As") == 0){return 3;}
			else if (s.compare("Co") == 0){return 1;}
			else if (s.compare("Cr") == 0){return 2;}
			else if (s.compare("Te") == 0){return 2;}
			else if (s.compare("Fe") == 0){return 0;}
			else if (s.compare("Pb") == 0){return 4;}
			else if (s.compare("Mn") == 0){return 1;}
			else if (s.compare("Hg") == 0){return 4;}
			else if (s.compare("Mo") == 0){return 2;}
			else if (s.compare("Ni") == 0){return 2;}
			else if (s.compare("Se") == 0){return 2;}
			else if (s.compare("Ti") == 0){return 3;}
			else if (s.compare("Zn") == 0){return 2;}
   			else if (s.compare("Si") == 0){return 4;}
			else if (s.compare("?") == 0){return 1;}
			else if (s.compare("Mg") == 0){return 2;}
			else if (s.compare("Li") == 0){return 1;}
			else if (s.compare("Al") == 0){return 3;}
			else if (s.compare("Zr") == 0){return 2;}
			 else if (s.compare("Bi") == 0){return 3;}
			else if (s.compare("Pd") == 0){return 2;}
			else if (s.compare("Pt") == 0){return 2;}
			else if (s.compare("Ru") == 0){return 3;}
			else if (s.compare("Rh") == 0){return 2;}			
			else if (s.compare("Ga") == 0){return 3;}
			else if (s.compare("Ge") == 0){return 4;}	
			else { 	
				return 3;			
				//cout << "Unknown symbol " << s << "\n";
				//exit(0);
			}		
		}


		int Molecule::getValance2(std::string s) {
			if (s.compare("S") == 0){return 6;}
//			else if (s.compare("I") == 0){return 7;}
			else { 				
				cout << "Unknown symbol " << s << "\n";
				cout << atoms << " " << bonds << "\n";
				for (int i = 0; i < atoms; i++) 
					if (s.compare(atomchar[i]) == 0) 
						cout << x[i] << " " << y[i] << "\n";
				exit(0);
			}		
		}
		
		
		void Molecule::setAtomMap(int atom, std::string neighbourhood) {
			//@agrume
			cout << "neighbourhood=" << neighbourhood << endl;
			int sumchars = 0;
			for (int t = 0; t < neighbourhood.length(); t++) {
				sumchars += (int)neighbourhood[t];					
			}
			
			atomMap[atom] = sumchars;
		}
		int Molecule::getBondType(int i,int j){
			return b_type[i][j];
		}
