#ifndef Node_h
#define Node_h 1

#include "General.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

class Node {
	private:
		int node_id;		
		Float x, y, z;

	public:
		int numH;
		int symbol, iso, charge, radical, stereo;
		int valance;
		int aroma;

		Node(int the_node_id, int the_symbol, int the_iso, int the_charge, int the_radical, int the_stereo, int the_valance, Float the_x, Float the_y, Float the_z, int the_aroma) :
		      node_id(the_node_id), symbol(the_symbol), iso(the_iso), charge(the_charge), radical(the_radical), stereo(the_stereo), valance(the_valance), x(the_x), y(the_y), z(the_z), aroma(the_aroma) 
		{	
		}

		void print() {
			cout << symbol << " " << iso << " " << charge << " " << radical << " " << stereo << " " << valance << " " << numH << "\n";
		};
		int getNodeID() {return node_id; };
		int getSymbol() {return symbol; };
		int getIso() {return iso; };
		int getCharge() {return charge; };
		int getRadical() {return radical; };
		int getStereo() {return stereo; };
		Float getX() {return x; };
		Float getY() {return y; };
		Float getZ() {return z; };

		void setHydrogen(int nh) { numH = nh; };
		~Node(){
                	//@debug
			//cout << "calling node destructor" << endl;
			////////
		}
};

#endif
