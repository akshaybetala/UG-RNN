#ifndef Molecule_h
#define Molecule_h 1

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "General.h"
#include <stack>
#include <queue>
#include <string>
#include "Graph.h"
#include "Node.h"

class Molecule {
	public:		
	
		int* atomMap;
		Node** all_atoms;
		int numrings;
		int** rings;
		int* ringsize;
		int* ra_num;	/* number of rings atom is a member of */
		int* ra_size;	/* size of ring atom is a mebember of */
		int* hyb;	/* hybridization */

		int* aroma; /* aromaticity flag */
		int* anti; /* anticlockwise chiral center flag*/
		int* hydrogenBonds;	
		int chiral;	/* global chiral */

		/* per atom */
		Float* x;
		Float* y;
		Float* z;	/* coordinates */
		int* u;	/* atom symbol */
		int* c;	/* atom charge */
		Float* mol2Charges;
		int* r; 	/* radical info, rare */
		int* iso;	/* isotope info, rare */
		int* chr;	/* chiral info */
		int* stereo;	/* stereo atom */
		int* valancy;
		int* hs;
		std::string* atomchar;
		std::string* atomletter;

		/* per bond */
	protected:	
		int** b_type;		/* bond connections and types */
	public:	
		int** b_stereo;	/* bond connections and stereo */

		int foo_int, i;
		std::string foo_string;

		int atoms;	/* number of atoms */
		int bonds;	/* number of bonds */
		
		Graph* g;	/* the only thing we need, the molecule as a graph with varoius operations */
		Float target;
		Float trueTarget;
		Float logP;
		Float X1;
		Float X2;
		Float mass;
		Float IC2;
		Float RB;
		Float* weights;

		std::string formula;

		Molecule(istream &is);
		~Molecule();
		//agrume
		int getBondType(int i,int j);
		int** getBondType(){return b_type;}
		/////////////////////////////
		void getMol2Info(istream &is);
		void loadweights();
		void writeformula();

	  	static int translateU(std::string s);
		int translateMol2U(std::string s);
		int getValance(std::string s);
		int getValance2(std::string s);
		void getCycles(int* u, int* k);

		void setAtomMap(int u, std::string s);

		void lightenMolecule();
};

#endif // Molecule_h
