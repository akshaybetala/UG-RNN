#include "Graph.h"
#include <stack>
#include <queue>
#include <string>

using namespace std;
//@agrume the_adjStereo=b_stereo, the_adjType=b_type
//@agrume the_total_vert=atoms_number, the_total_edges=bonds_number 
Graph::Graph(Node** the_vert, int** the_adjType, int** the_adjStereo, int the_total_vert, int the_total_edges) : 
		vert(the_vert),adjStereo(the_adjStereo), total_vert(the_total_vert), total_edges(the_total_edges)
{
	adjType = new int*[total_vert];
	for(int a=0;a<total_vert;a++){
		adjType[a] = new int[total_vert]; 
		for(int b=0;b<total_vert;b++){
			adjType[a][b]=the_adjType[a][b];
		}
	}
	

}

void Graph::orderGraph() {
	pairSP = new int*[total_vert];
	for (int i = 0; i < total_vert; i++) {
		pairSP[i] = new int[total_vert];
		memset(pairSP[i], 0, sizeof(int)*total_vert);
	}

	order = new int*[total_vert];
	for (int i = 0; i < total_vert; i++) {
		order[i] = new int[total_vert];
	}
	//cout << "Dijkstra \n " << flush;
	/* create the shortest path matrix using Dijkstra N times
	   Then order graph towards shortest path
	   Also direct equal paths with most info (e.g. Not Carbon atoms): Not done yet  */
	for (int i = 0; i < total_vert; i++) {
		dijkstra(i);
	}	
	//cout << "Creating topologiacal order \n" << flush;
	for (int i = 0; i < total_vert; i++) {
		orderGraph(i);
	}

	//cout << "Creating SCC \n" << flush;
	getScc();
	setMainScc();

	for (int a = 0; a < total_vert; a++) {		
		delete[] pairSP[a];		
	}
	delete[] pairSP;


	/* print adj matrix 
	for (int i = 0; i < total_vert; i++) {
		for (int j = 0; j < total_vert; j++) {
			cout << adjType[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";
	*/
}

void Graph::loadNeighbours() {
	//cout << "Finding children and parents for each order \n" << flush;
	children = new int**[total_vert];
	parents = new int**[total_vert];
	num_children = new int*[total_vert];
	num_parents = new int*[total_vert];
	for (int i = 0; i < total_vert; i++) {
		children[i] = new int*[total_vert];
		parents[i] = new int*[total_vert];
		for (int j = 0; j < total_vert; j++) {
			children[i][j] = new int[MAXC];
			parents[i][j] = new int[MAXC];
			for (int c = 0; c < MAXC; c++) {
				children[i][j][c] = -1;
				parents[i][j][c] = -1;
			}
		}
		num_children[i] = new int[total_vert];
		memset(num_children[i], 0, sizeof(int)*total_vert);
		num_parents[i] = new int[total_vert];
		memset(num_parents[i], 0, sizeof(int)*total_vert);
	}
	
	//cout << "Setting neighbours \n" << flush;
	for (int i = 0; i < total_vert; i++) {
		//cout << "neighbours for " << i << "\n" << flush;
		setNeighbours(i);
	}
}

void Graph::unloadGraph() {
	//cout << "Deleting Graph (actually " << scc  << " graphs) " << "\n" << flush;
	for (int a = 0; a < total_vert; a++) {
		for (int j = 0; j < total_vert; j++) {
			delete[] parents[a][j];
			delete[] children[a][j];
		}
		delete[] parents[a];
		delete[] children[a];
		delete[] num_parents[a];
		delete[] num_children[a];
		delete[] order[a];
	}
	delete[] parents;
	delete[] children;
	delete[] num_parents;
	delete[] num_children;
	delete[] order;
	for(int s = 0; s < scc; s++) {		
		for (int v = 0; v < scclength[s]; v++) {
			delete[] scc_order[s][v];
		}
		delete[] scc_order[s];
	}
	delete[] scc_order;
	delete[] scclength;
}

void Graph::lightenGraph() {
	//cout << "Deleting Graph (actually " << scc  << " graphs) " << "\n" << flush;
	for (int a = 0; a < total_vert; a++) {
		for (int j = 0; j < total_vert; j++) {
			delete[] parents[a][j];
		}
		delete[] parents[a];
		delete[] num_parents[a];
		delete[] order[a];
	}
	delete[] parents;
	delete[] num_parents;
	delete[] order;
	for(int s = 0; s < scc; s++) {		
		for (int v = 0; v < scclength[s]; v++) {
			delete[] scc_order[s][v];
		}
		delete[] scc_order[s];
	}
	delete[] scc_order;
	delete[] scclength;
}



Graph::~Graph() {
	for (int a = 0; a < total_vert; a++) {
		delete[] adjType[a];
	}
	
	delete[] adjType;

}


/* Sets the children and the parents of each node given a supersource */
/* O(who cares) it only happens at begining of training */
void Graph::setNeighbours(int root) {

	int curnode;
	int curnode_child;
	int curindex;
	for (int i = 0; i < total_vert; i++) {
		curnode = i;
		/* Where in the order is curnode */
		curindex = -1;
		for (int n = 0; n < total_vert; n++) {
			if (curnode == order[root][n]) {
				curindex = n;
				break;
			}
		}
		if (curindex >= 0) {	/* otherwise the current node is not part of the graph with this root */

			for (int j = 0; j < total_vert; j++) {
				if (adjType[i][j] > 0) {	/* a neighbour */
					curnode_child = j;
					/* Now see if curnode is left (child) or right (parent) in ordered graph */
					for (int n = curindex; n >=0; n--) {
						if (curnode_child==order[root][n]) {
							children[root][curnode][num_children[root][curnode]] = curnode_child;
							num_children[root][curnode]++;
							//cout << "root " << root << "node now " << curnode << " Cnumsofar " << num_children[root][curnode] << "\n" << flush;
							break;
						}
					}
					for (int n = curindex; n < total_vert; n++) {
						if (curnode_child==order[root][n]) {
							parents[root][curnode][num_parents[root][curnode]] = curnode_child;
							num_parents[root][curnode]++;
							//cout << "root " << root << "node now " << curnode << " Pnumsofar " << num_parents[root][curnode] << "\n" << flush;
							break;
						}
					}
					if (num_children[root][curnode]+num_parents[root][curnode] > MAXC) {
						cout << "WARNING: NODE " << curnode << " has " << (num_children[root][curnode]+num_parents[root][curnode]) << " neighbours (calculated so far), max allowed is " << MAXC << "\n" << flush;
						//num children
						cout << "num_children" << num_children[root][curnode] << endl;
						exit(0);
					}
				}
			}
		}
	}
/*
	if (scc == 2) {
	for (int i = 0; i < total_vert; i++) {
		cout << "Direction towards " << root << " children[" << i << "]=" << num_children[root][i] << "\n";
		cout << "                    parent[" << i << "]=" << num_parents[root][i] << "\n";
		for (int c = 0; c < num_children[root][i]; c++) 
			cout << children[root][i][c] << " ";
		cout << "\n";
	}
	}
*/

}
/* Get main graph from all the scc's. 
   Simplest: Just get biggest graph */
void Graph::setMainScc() {
	int max;
	if (scc !=0) {
		max = 0;
		for (int s = 0; s < scc; s++) {
			if (scclength[s] > max) {
				max = scclength[s];
				mainscc = s;
			}
		}
	}
	else {mainscc = -1;}
}

/* Extract all the Graph's from the group of 
graphs. The graphs are separated using the node order matrix.
I think the groups are called strongly connected components (fancy ;->). 
*/
void Graph::getScc() {
	scc = 0;	/* global*/

	int jumpto = 0;
	int* tmp = new int[2000];
	int* done = new int[total_vert];
	memset(done, 0, sizeof(int)*total_vert);	
	memset(tmp, 0, sizeof(int)*2000);	
	int* vertex = new int[2000];
	vertex[0] = 0;
	int alldone = 0;
	while (!alldone) {		
		for (int v2 = 0; v2 < total_vert; v2++) {
			if (order[vertex[scc]][v2] >= 0) {
				tmp[scc]++;
				done[order[vertex[scc]][v2]]=1;
			}
		}
		for (int v2 = 0; v2 < total_vert; v2++) {
			if (!done[v2]) {
				scc++;
				vertex[scc] = v2;
				alldone = 0;
				break;
			}
			else {
				alldone = 1;
			}
		}		
	}
	scc++;
	scc_order = new int**[scc];
	scclength = new int[scc];		/*global */
	for (int s = 0; s < scc; s++) {
		scc_order[s] = new int*[tmp[s]];
		for (int v = 0; v < tmp[s]; v++) {
			scc_order[s][v] = new int[tmp[s]];		
		}
		scclength[s] = tmp[s];
	}
	delete[] tmp;
	
	int index = 0;
	for(int s = 0; s < scc; s++) {		
		int* todo = new int[scclength[s]];
		memset(todo, 0, sizeof(int)*scclength[s]);	
		index = 0;
		for (int v = 0; v < total_vert; v++) {
			if (order[vertex[s]][v] >= 0) {
				todo[index] = order[vertex[s]][v];
				index++;
			}
		}
		for (int v = 0; v < scclength[s]; v++) {
			for (int v2 = 0; v2 < scclength[s]; v2++) {
				scc_order[s][v][v2] = order[todo[v]][v2];			
			}
		}
		delete[] todo;
	}

	delete[] vertex;
	delete[] done;
}

void Graph::getCycles() {

	/* Calculate all shortest cycles. These can be reduced to one node using reduceGraph() */
	cyclepoints = new int*[total_vert];
	for (int v = 0; v < total_vert; v++) {
		cyclepoints[v] = new int[total_vert];
		memset(cyclepoints[v], 0, sizeof(int)*total_vert);
	}
	numcycles = new int[scc];
	memset(numcycles, 0, sizeof(int)*scc);
	cycles = new int**[scc];
	for (int s = 0; s < scc; s++) {
		for (int v = 0; v < total_vert; v++)
			//?
			memset(cyclepoints[v], 0, sizeof(int)*total_vert);

		cout << "Calculating Back Edges.. \n";
		dfs(s);

		cout << "Number of cycles " << numcycles[s] << "\n";
		numcycles[s] = 20;	
		cycles[s] = new int*[numcycles[s]];
		for (int c = 0; c < numcycles[s]; c++) {
			cycles[s][c] = new int[total_vert];
			for (int v = 0; v < total_vert; v++) {
				cycles[s][c][v] = -1;
			}
		}
		
		int c = 0;
		for (int i = 0 ; i < total_vert; i++) {
			for (int j = 0 ; j < total_vert; j++) {
				if (cyclepoints[i][j] > 0 && c < numcycles[s]) {
					c = shortestCycle(i,j,c,s);
				}	
			}
		}
		cout << "Cycles found::\n";
		for (int c = 0; c < numcycles[s]; c++) {
			for (int i = 0 ; i < total_vert; i++) {
				cout << cycles[s][c][i] << " ";
			}
			cout << "\n";
		}
	}

	for (int a = 0; a < total_vert; a++) {
		delete[] cyclepoints[a];
	}
	delete[] cyclepoints;
}

int Graph::shortestCycle(int start, int end, int cy, int sc) {
	int min, next;
	int soFar = 0;
	for (int dos = 0; dos < sc; dos++) 
		soFar += scclength[dos];

	int* d = new int[total_vert];
	int* s = new int[total_vert];

	cout << "Finding shortest cycles between " << start << " " << end << "\n";
	for (int v = 0; v < total_vert; v++) {
		d[v] = 10000;
		s[v] = 0;
	}
	int child = -1;
	for (int v = scclength[sc]-2; v >= scclength[sc]-1-getNumChildren(start); v--) {
		child = scc_order[sc][start-soFar][v];
		d[child] = 1;
	}
	d[end] = cyclepoints[start][end]+1;
	d[start] = 0;
	s[start] = 1;

	/* Assume that there at maximum two cycles for one edge. Is This TRUE?? */	
	int* u = new int[scclength[sc]];
	int* k = new int[scclength[sc]];
	int i;
	int n;
	int foundend;
	int distinctcycle;
	int whichcycle;	
	if (cy==0)
		distinctcycle = 1;
	else
		distinctcycle = 0;
	for (int docycles = 0; docycles < 2; docycles++) {
		for (int j = 0; j < scclength[sc]; j++) 
			u[j] = -1;
		for (int j = 0; j < scclength[sc]; j++) 
			k[j] = -1;
		s[end] = 0;		
		i = 0;
		n = 0;
		foundend = 0;
		while (!foundend && i < scclength[sc]) {
			next = 10000;
			min = 10000;
			for (int j = 0; j < total_vert; j++) {
				if (!s[j]) {		
					if (d[j] < min) {
						min = d[j];
						next = j;
					}
				}
			}
			if (next == 10000) {
				//cout << "Nowhere to go \n";
			}
			else {				
				s[next] = 1;
				d[next] = min;
				for (int v = scclength[sc]-2; (v >= scclength[sc]-1-getNumChildren(next)) && !foundend; v--) {
					child = scc_order[sc][next-soFar][v];
					if (d[next] + 1 < d[child]) {
						d[child] = d[next] + 1;	
						u[n] = next;
						k[n] = child;	
//						cout << u[n] << " " << k[n] << "\n";
						n++;
						if (child==end && n>1)
							foundend = 1;
						//cout << next << " " << v << "\n";
					}
				}
			}
			i++;
		}
	
		int cur, min, min_index, tmp, foundcycle;
		if (!foundend) {
			cout << "Cannot find cycle YET\n";
		}
		else {	
			/* where to place the cycle */
			for (whichcycle = 0; whichcycle < numcycles[sc]; whichcycle++) {
				if (cycles[sc][whichcycle][0] == -1) { break; }
			}

			/* now lets load the cycle from start to end and block the cycle for the future */
			for (int j = 0; j < total_vert; j++) 	
				s[j]=0;			
			n--;
			int moveto;
			int looplength;
			int cl = 0;
			cycles[sc][whichcycle][0] = start;
			cycles[sc][whichcycle][1] = k[n];
			cycles[sc][whichcycle][2] = u[n];
			cl = 3;
			while (n >=0) {
				moveto = u[n];
				while (moveto != k[n] && n>=0) {			
					n--;
				}	
				if (n >= 0) { 
					looplength++;
					cycles[sc][whichcycle][cl] = u[n];
					cl++;
				}
			}
//			cycles[sc][cy][cl] = start;
			s[cycles[sc][whichcycle][cl/2]] = 1;	/* cl/2 is the furthest vertex from backedge (start,end), this is done to block path to already seen cycle. MAY NOT WORK ALWAYS*/
			cout << "blocking vertex " << cycles[sc][whichcycle][cl/2] << "\n";
			
			/* restart distances */
			for (int v = 0; v < total_vert; v++) {
				d[v] = 10000;		
			}
			int child = -1;
			for (int v = scclength[sc]-2; v >= scclength[sc]-1-getNumChildren(start); v--) {
				child = scc_order[sc][start-soFar][v];
				d[child] = 1;
			}
			d[end] = cyclepoints[start][end]+1;
			d[start] = 0;
			s[start] = 1;
			cout << "Curently blocked edges \n";
			for (int v = 0; v < total_vert; v++) {
				cout << s[v] << " ";
			}
			cout << "\n";
			cout << "Before Sort: The cycle found and the position placed " << whichcycle << "\n";
			for (int v = 0; v < total_vert; v++) {
				cout << cycles[sc][whichcycle][v] << " ";
			}
			cout << "\n";
			/* sort the cycle to make it easy to see if already there later */
			for (int v = 0; v < cl; v++) {
				min = 1000;
				min_index=-1;
				for (int v2 = v; v2 < cl; v2++) {
					if (cycles[sc][whichcycle][v2]>=0)
						cur = cycles[sc][whichcycle][v2];
					if (cur < min) {
						min = cur;
						min_index = v2;
					}
				}		
				tmp = cycles[sc][whichcycle][v];
				cycles[sc][whichcycle][v] = min;
				cycles[sc][whichcycle][min_index] = tmp;
			}
			cout << "After Sort: The cycle found and the position placed " << whichcycle << "\n";
			for (int v = 0; v < total_vert; v++) {
				cout << cycles[sc][whichcycle][v] << " ";
			}
			cout << "\n";

			/* Have we already found the cycle? */
			foundcycle = 1;			
			for (int c = 0; (c < numcycles[sc]); c++) {
				if (c!=whichcycle && cycles[sc][c][0]!=-1) {
					foundcycle = 1;
					for (int v = 0; v < total_vert; v++) {					
						if (cycles[sc][whichcycle][v] != cycles[sc][c][v]) {
							foundcycle = 0;
						}
					}
					if(foundcycle) {
						cout << "Found cycle already. Deleting..... \n";
						for (int v = 0; v < total_vert; v++) {				
							cycles[sc][whichcycle][v] = -1;
						}
					}
				}
			}
		}
	}
	delete[] d;
	delete[] s;
	delete[] u;
	delete[] k;	
	cout << whichcycle << " cycles found SOFAR\n\n";
	return whichcycle;
}

/* Shortest path from node root. This will also 
   determine the number of parents and children of each node
   at each supersource */
void Graph::dijkstra(int root) {	
	int min, next;
	int* d = new int[total_vert];
	int* s = new int[total_vert];
	for (int v = 0; v < total_vert; v++) {
		if (adjType[root][v] > 0) {
			d[v] = 1;		/* atomic dist?? */
		}
		else {
			d[v] = 10000;
		}
		s[v] = 0;
	}
	d[root] = 0;
	s[root] = 1;

	for (int n = 1; n < total_vert; n++) {
		next = 10000;
		min = 10000;

		/* get min vert, u.  */
		for (int j = 0; j < total_vert; j++) {
			if (!s[j]) {		
				if (d[j] < min) {
					min = d[j];
					next = j;
				}
			}
		}
		if (next == 10000) {
			//cout << "Distinct graph \n";
		}
		else {
			s[next] = 1;
			d[next] = min;
			for (int v = 0; v < total_vert; v++) {
				if (adjType[next][v] > 0) {
					if (d[next] + 1 < d[v]) {
						d[v] = d[next] + 1;					
					}
				}
			}
		}
	}
	for (int v = 0; v < total_vert; v++) {
		if (d[v] != 10000) {
			pairSP[root][v] = d[v];
			pairSP[v][root] = d[v];
		}
		else {
			pairSP[root][v] = -1;
			pairSP[v][root] = -1;
		}
	}
/*
	for (int v = 0; v < total_vert; v++) 
		cout << pairSP[root][v] << " ";
	cout << "\n";
*/
	delete[] d;
	delete[] s;
}

/* This is vital for processing the graph. It shows the visit order of the graph.
   From the furtest vertex to the vertex currently being examined.
*/
void Graph::orderGraph(int root) {
	int* visit = new int[total_vert];

	int max = 0;
	int max_node = -1;

	for (int v = 0; v < total_vert; v++) {
		if(pairSP[root][v] == 10000) 
			visit[v] = 1;	/* disjoint graph. Lets not consider for the moment */
		else
			visit[v] = 0;
	}
	/* O(n^2) sort i dont mind only happens at beigning of training */
	for (int v = 0; v < total_vert; v++) {
		max = 0;
		max_node = -1;
		for (int i = 0; i < total_vert; i++) {
			if (!visit[i]) {
				if (pairSP[root][i] >= max) {
					max = pairSP[root][i];
					max_node = i;
				}
			}
		}
		if (max_node >= 0) {
			order[root][v] = max_node;
			visit[max_node] = 1;
		}
		else {
			order[root][v] = -1;	/* -1 for disjoint graph */
		}
	}

	/*
	cout << "root " << root  << "\n" << flush;
	for (int v = 0; v < total_vert; v++) 
		cout << order[root][v] << " ";
	cout << "\n";
	*/
	delete[] visit;
}



int Graph::isCycle() {
		
	return 0;
}

void Graph::dfs(int sc) {
	/* get the number of nodes processed so far */
	int soFar = 0;
	for (int dos = 0; dos < sc; dos++) 
		soFar += scclength[dos];

	std::stack<int> s;
	int* dfsorder = new int[total_vert];

	int* visit = new int[total_vert];
	memset(visit, 0, sizeof(int)*total_vert);


	int curV;
	int i = 0;
	s.push(scc_order[sc][0][scclength[sc]-1]);
	int* camefrom = new int[total_vert];
	int* stackelements = new int[1000];
	int alreadythere = 0;
	int j = 1;
	while(s.size()!=0 && i < scclength[sc]) {	/* fix this to main molecule */

		curV = s.top();	
		s.pop();
		j--;
		visit[curV]=1;
		dfsorder[i] = curV;
//		cout << curV << "\n";
		int child = -1;
		alreadythere = 0;
		for (int v = scclength[sc]-2; v >=scclength[sc]-1-getNumChildren(curV); v--) {
			child = scc_order[sc][curV-soFar][v];			
			if (!visit[child]) {				
				for (int n = 0; n < j; n++) {
					if (stackelements[n]==child)
						alreadythere=1;
				}
				if (!alreadythere) {
					s.push(child);
					stackelements[j] = child;
					j++;
					camefrom[child] = curV;
				}
			}
			else {
				if (camefrom[curV] != child) {
					if (cyclepoints[child][curV]!=10000) {
						cout << curV << " " << child << " " << camefrom[curV] << "\n";						
						cyclepoints[curV][child] = 10000;
//						cyclepoints[child][curV] = 10000;
						numcycles[sc]++;
					}						
				}
			}
		}
		i++;		
	}
				
	delete[] stackelements;
	delete[] dfsorder;
	delete[] visit;
	delete[] camefrom;
}


Node* Graph::getNode (int i) {
	return vert[i];
}


