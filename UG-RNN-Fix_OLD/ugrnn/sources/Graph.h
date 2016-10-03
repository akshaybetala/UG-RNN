#ifndef Graph_h
#define Graph_h 1

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "Node.h"
#include "General.h"

#define MAXC 40
#define MAX_CYCLES 10

class Graph {

	private:		
		int** adjStereo;

		int** pairSP;
		int** order;			
		int* numcycles;
		int** cyclepoints;
		int*** cycles;

	public:
		int** adjType;
		Node** vert;					/* following variables, used throughout training */
		int total_vert;	
		int total_edges;
		int scc;

		int** num_children;	/* Depends on the supersource, thats why 2D*/
		int** num_parents;	
		int*** children;
		int*** parents;

		int* scclength;
		int*** scc_order; /* visit order of all strongly connected components */
		int mainscc;

		Graph(Node** the_vert, int** the_adjType, int** the_adjStereo, int the_total_vert, int the_total_edges);
		~Graph();

		int getTotalVert() { return total_vert; }
		int getTotalEdges() { return total_edges; }
		int getNumChildren(int v) { return num_children[0][v]; }	/* TO DO */

		void orderGraph(int root);	/* direct graph to root */
		void dijkstra(int root);
		int shortestCycle(int root, int end, int cycle, int sc);
		void getCycles();
		void getScc();
		Node* getNode(int i);

		int isCycle();

		void reduceGraph(); 	/* reduce graph by choosing larger atoms (e.g. let a ring be an atom) */		
		void setMainScc();
		void setNeighbours(int root);
		void getBonds();

		void dfs();
		void dfs(int s);
		void dfs(int start, int end);
		void bfs();

		void orderGraph();
		void unloadGraph();
		void lightenGraph();
		void loadNeighbours();
};

#endif
