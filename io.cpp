//#include "gsjoin.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

int gsjoin::ReadParameters(char* sim_func, char* threshold)
{
	if (sim_func[0] == 's' || sim_func[0] == 'S') ss = 1;
	else return 1;
	gsjoin::threshold = atof(threshold);
	cerr << "Similarity Function: ";
	if (ss == 1) cerr << "Sequence Similarity" << endl;
	cerr << "Threshold: " << gsjoin::threshold << endl;
	return 0;
}


void gsjoin::WriteResults()
{
	int i;
	//ofstream outfile("Results0.8.txt");
	vector<result>::iterator it;
	for (it = joinResult.begin(); it != joinResult.end(); ++it)
		cout << it->xid << " " << it->yid << " " << setiosflags(ios::fixed) << setprecision(3) << it->sim << endl;
		//outfile << it->xid << " " << it->yid << " " << setiosflags(ios::fixed) << setprecision(3) << it->sim << endl;
		
}

//read graph data for VR similarity measure
void gsjoin::readData(char *filename)
{
	int c, edges;
	FILE *fp;
	fp=fopen(filename,"r");
	fscanf(fp,"%d",&c);
	for(int j=0;j<c;j++)
	{
		int gid,edgeCount;
		fscanf(fp,"%d",&gid);
		fscanf(fp,"%d",&edgeCount);
		cout << gid << endl;
		vector< pair<int,int> > graph;
		for(int i = 0; i < edgeCount; i++)
		{
			int from,to;
			fscanf(fp,"%d %d",&from,&to);
			cout<<from<<" "<<to<<endl;
			graph.push_back(make_pair(from,to));
		}
		idGraph[gid]=graph;
	}
	fclose(fp);
}

//reading a graph in new data set format  (----------- as delimiter)
void gsjoin::readDataNew(char *filename)
{
	//declare an input file stream
	ifstream inFile;
	string x;
	//open the file stream
	inFile.open(filename);
	int gid = 0; //used to index the graphs
	//creating array of graphs from input file
	while (!inFile.eof()) {
	  	//read  a line from the file every time
		getline(inFile,x);
		//if the line is containing an edge of a graph
		if(isdigit(x[0]))
		{
			stringstream stream(x);
			int a,b;
			stream >> a >>  b ; //get the co-ordinates for the edge from string stream
			//add the edge to the edgelist of the graph
			idGraph[gid].push_back(make_pair(a,b)); 
		}
		//if the line is containing delimiter (--------)
		else
		{
			if(x.compare("----------") == 0)
			{	
				gid++;
				//create a new empty vector of pairs to insert edges
				vector< pair<int,int> > edgeList;
				idGraph[gid] = edgeList;
			}
		}
	}
}
