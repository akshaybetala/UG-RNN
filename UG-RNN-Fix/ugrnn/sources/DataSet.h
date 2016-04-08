#include <iostream>
#include "General.h"
#include "Molecule.h"
#include <queue>
#include <string>
#include <sstream>

using namespace std;

class DataSet {

	public:
		int length;
		Molecule** mole;
		
		int max_ringnumber;
		int max_ringsize;
		int noLogP;
		Float pos;
		Float neg;
		Float skew;
		Float classes;
		Float max;
		Float min;
		Float myMax;
		Float myMin;
		Float maxLogP;
		Float minLogP;
		vector<Float> trueTargets;

		DataSet() {}

		DataSet(istream& is, int classes) {
			//sym = new int[100];
			//memset(sym, 0, sizeof(int)*100);
			
			is >> length;
			cout << "agrume: number of molecules in dataset = " << length << endl;
		        //@agrume number of molecules in dataset
			mole = new Molecule*[length];
			max_ringnumber = 0;
			max_ringsize = 0;
			pos = 0;
			neg = 0;
			trueTargets = vector<Float>(length,0);
			cout << "LOADING MOLECULES..." << endl;
			for (int i = 0; i < length; i++) {
			//for (int i = 0; i < length; 1) {
				//@debug
				cout << "loading molecule " << i << flush << endl;
				////////
				mole[i] = new Molecule(is);
				mole[i]->trueTarget = mole[i]->target;
				mole[i]->lightenMolecule();	
				//debug
				//cout << mole[i]->target << endl;
				//cout << mole[i]->logP << endl;
				///////			
				//trueTargets[i] = mole[i]->target;
				//exit(0);
			}
			cout << "...DONE" << endl;
		}

		
		void normalizeTargets(Molecule** molecule,int arrayLength){
			Float sum = 0.0;
			Float avg = 0.0;
		 	max = molecule[0]->target;
		 	min = molecule[0]->target;
			maxLogP = molecule[0]->logP;
			minLogP = molecule[0]->logP;
			myMin = -0.5;
		 	myMax = 0.5;			
			Float myRange = myMax-myMin;		 	
			Float std = 0;
			noLogP = 1;
			cout << "calcuating mean and std \n" << flush;
			for (int i = 0; i < arrayLength; i++) {
				cout << "." << flush;
				sum += molecule[i]->target;
				if (molecule[i]->target > max)
					max = molecule[i]->target;
				if (molecule[i]->target < min)
					min = molecule[i]->target;
				if (molecule[i]->logP > maxLogP)
					maxLogP = molecule[i]->logP; 
				if (molecule[i]->logP < minLogP)
					minLogP = molecule[i]->logP;
				if (molecule[i]->logP != 0.0)			
					noLogP = 0;			
			}
			cout << "\n";
			avg = sum/(double)arrayLength;
			std = 0;
			for (int i = 0; i < arrayLength; i++) {
				std += (molecule[i]->target-avg)*(molecule[i]->target-avg);
			}
			std = std/(double)arrayLength;
			std = sqrt(std);
			for (int i = 0; i < arrayLength; i++) {
			//molecule[i]->target = mole[i]->target - avg;
				Float target = molecule[i]->target;
				cout << molecule[i]->target << " ";
				molecule[i]->target = molecule[i]->target - min;
				//molecule[i]->target = mole[i]->target / (Float)std;
				molecule[i]->target = (Float)myMin+(Float)myRange*(molecule[i]->target / (Float)(max-min));
				//cout << (Float)(molecule[i]->target-myMin)*(Float)(max-min)/(Float)(myMax-myMin)+(Float)min << endl;
				cout << molecule[i]->target << endl;
				//debug
				//cout << molecule[i]->target << "\n";
				///////
				//normalize logP
				//cout << "maxLogP=" << maxLogP << " " << "minLogP=" << minLogP << endl;
				if (noLogP==0)				
					molecule[i]->logP = (Float)myMin+(Float)myRange*(molecule[i]->logP / (Float)(maxLogP-minLogP));
				////////////////	
			}
			//cout << "Calculated mean : " << avg << " calculated std : " << std << "\n";
			sum = 0.0;
			avg = 0;
			std = 0;
			for (int i = 0; i < arrayLength; i++) {
				sum += molecule[i]->target;
			}
			avg = sum/(double)arrayLength;
			std = 0;
			for (int i = 0; i < arrayLength; i++) {
				std += (molecule[i]->target-avg)*(molecule[i]->target-avg);
			}
			std = std/(double)arrayLength;
			std = sqrt(std);
			cout << "new_mean=" << avg << " new_std=" << std << "\n";
			cout << "min=" << min << " max=" << max << endl;
			cout << "noLogP=" << noLogP << endl;				
		}
		
		
		//edited by agrume	
		void normalizeTargets(Molecule** molecule,int arrayLength,Float ma,Float mi){
			Float sum = 0.0;
			Float avg = 0.0;
			myMin = -0.5;
			myMax = 0.5;
			Float myRange = myMax-myMin;
			Float std = 0;
			cout << "calcuating mean and std \n" << flush;
        	    	for (int i = 0; i < arrayLength; i++) {											
				cout << "." << flush;
				sum += molecule[i]->target;
                        }
			cout << "\n";
			avg = sum/(double)arrayLength;                                
			std = 0;
			for (int i = 0; i < arrayLength; i++) {
				std += (molecule[i]->target-avg)*(molecule[i]->target-avg);
			}
			std = std/(double)arrayLength;
			std = sqrt(std);
			for (int i = 0; i < arrayLength; i++) {
				//molecule[i]->target = mole[i]->target - avg;
				Float target = molecule[i]->target;
				cout << molecule[i]->target << " ";
				molecule[i]->target = molecule[i]->target - mi;
				molecule[i]->trueTarget = target;
				//molecule[i]->target = mole[i]->target / (Float)std;
				molecule[i]->target = (Float)myMin+(Float)myRange*(molecule[i]->target / (Float)(ma-mi));
				//cout << (Float)(molecule[i]->target-myMin)*(Float)(ma-mi)/(Float)(myMax-myMin)+(Float)mi << endl;
				cout << molecule[i]->target << endl;
				//@debug
				//cout << molecule[i]->target << "\n";
				////////
				if(noLogP==0)	
					molecule[i]->logP = (Float)myMin+(Float)myRange*(molecule[i]->logP / (Float)(maxLogP-minLogP));	
			}
			//@debug
			//cout << "Calculated mean : " << avg << " calculated std : " << std << "\n";
			sum = 0.0;
			avg = 0;
			std = 0;
			for (int i = 0; i < arrayLength; i++) {
				sum += molecule[i]->target;
			}
			avg = sum/(double)arrayLength;                               
			std = 0;
			for (int i = 0; i < arrayLength; i++) {
				std += (molecule[i]->target-avg)*(molecule[i]->target-avg);
			}
			std = std/(double)arrayLength;
			std = sqrt(std);
			cout << "new_mean=" << avg << " new_std=" << std << "\n";
			cout << "min=" << mi << " max=" << ma << endl;
			cout << "noLogP=" << noLogP << endl;		
		}
		
		void deNormalizeTargets(){
			for(int i=0;i<length;i++){
				mole[i]->target = trueTargets[i];
				//debug
				cout << mole[i]->target << endl;
				///////
			}	
		}	
};
