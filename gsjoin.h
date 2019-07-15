#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <bits/stdc++.h>
#include <map>

using namespace std;

struct result
{
	int xid, yid;
	double sim;
};


class gsjoin {
	bool ss;
	double threshold;
	int numHash;
	vector<result> joinResult;
	map< int,vector< pair<int,int> > > idGraph;	// GID and its edges
	map< int,vector< pair<int,double> > > idPageRank; //GID and a pair contians vertex and its PageRank
	map<int, vector< pair<int,int> > > sortedGraphs; //GID and a pair contains vertex and its ranking based its PageRank
	map<int, vector< pair<int,int> > > sortedGraphsOld; //used to build the sorted graphs
	map<int, set<int> > graphShingles;
	map<int, set<int> > invertedIndex ; // Shingle and A Set of Graph id'S containing that single
	map<int, set<int> > relatedGaphs; // GID and a set of Graph ID's which are related it


	set<int> unionShingles;
	int** hashCoeffs;
	public:
		gsjoin() {
			ss = 0;
			numHash = 5;
			unionShingles.clear();
			hashCoeffs = new int*[numHash];
			for(int i = 0; i < numHash; ++i)
   				 hashCoeffs[i] = new int[2];

		}
		~gsjoin()
		{
			for(int i = 0; i < numHash; ++i) {
    			delete [] hashCoeffs[i];
			}
			delete [] hashCoeffs;
		}
		int ReadParameters(char*, char*);
		void readData(char*);
		void readDataNew(char*);
		void getPageRanks();
		void rankingVertex();
		static bool sortbysec(const pair<int,double>&, const pair<int,double>&);
		void getShingles();
		void computeUnionShingles();
		void computeHashCoeffs();
		unsigned int computeHash(int, int, int);
		int getMin(int, int);
		vector<int> getBoolShingles(set<int>&);
		vector<unsigned int> getGraphMinHash(vector<int>&, int);
		void printGraphShingles(set<int>&);
		void printGraphShinglesBool(vector<int>&);
		void printGraphSignature(vector<unsigned int>&);
		void performJoinSequenceSimilarity();
		void WriteResults();

		void rankingVertexNew(); //to rank the vertex by considering adjacent vertices and pagerank
		void  getSortedSequence(int);//to get the sorted sequence for the given graph representation (set of edges)
		vector<int> getAdjList(int,int); //to get adjacency list of vertex for given graph 		
		void performJoinSequenceSimilarityNew(); //based on the inverted index 
		void constructInvertedIndex(); // to construct inverted index from given set of graphs
};
