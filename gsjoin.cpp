#include <math.h>
//#include <google/dense_hash_map>
#include <algorithm>
#include <bits/stdc++.h>
#include <map>
#include <vector>
#include <iterator>
#include <utility>
#include <cstdlib>
#include <climits>
#include "gsjoin.h"
//#include "index.cpp"
//#include "io.cpp"
using namespace std;


void gsjoin::getPageRanks()
{
	//cout<<"Dumping edges\n";
	map<int,vector< pair<int,int> > >::iterator it;
	for(it=idGraph.begin();it!=idGraph.end();it++)
	{
		//cout<<it->first<<endl;
		vector< pair<int,int> > graph=it->second;
		int gid=it->first;		
		FILE *fp;
		fp=fopen("dataGraph.txt","w");
		//dataGraph.txt contains edges of one graph at a time
		for(int i=0;i<graph.size();i++){
			//printf("%d %d\n",graph[i].first,graph[i].second);
			fprintf(fp,"%d %d\n",graph[i].first,graph[i].second);
		}
		fclose(fp);
		const int BUFSIZE = 1000;
    	char buf[ BUFSIZE ];
    	//calling pageRank algorithm for each graph
		FILE * f = popen( "./pg -d ' ' dataGraph.txt", "r" );
		while( fgets( buf, BUFSIZE,  f ) ) {
        		fprintf( stdout, "%s", buf  );
    		}
    		pclose(f);
    		//sleep(1);
    		fp=fopen("out.txt","r");
    		if(fp==NULL){
    			printf("File reading Error !!\n");
    			exit(0);
    		}
    		int vertex;
    		double pagerank;
    		//pageRankGraph contains vertex and its pageRank score
    		vector<pair<int,double> > pageRankGraph;
    		while(fscanf(fp,"%d %lf\n",&vertex,&pagerank)!=EOF){
    			pageRankGraph.push_back(make_pair(vertex,pagerank));
    		}
    		idPageRank[gid]=pageRankGraph;
    		fclose(fp);
    		//reading vertex and its #outlinks
    		fp = fopen("outlinks.txt","r");
    		if(fp==NULL){
    			printf("File reading Error !!\n");
    			exit(0);
    		}
}
	//printing all graphs :vertices and their pageranks
	/*for(map< int,vector< pair<int,double> > >::iterator it=idPageRank.begin();it!=idPageRank.end();it++ ){
		cout<<it->first<<endl;
		vector< pair<int,double> > graphPageRank = it->second;
		for(int i=0;i<graphPageRank.size();i++){
			cout<<graphPageRank[i].first<<" "<<graphPageRank[i].second<<endl;
		}
		cout<<endl;
	}
	//printing vertices and their outlinks for each graph
	for(map< int,vector< pair<int,int> > >::iterator it = graphOutlinks.begin();it != graphOutlinks.end(); it++ ){
		cout<<it->first<<endl;
		vector< pair<int,int> > vertexOutlinks = it->second;
		for(int i=0;i<vertexOutlinks.size();i++){
			cout<<vertexOutlinks[i].first<<" "<<vertexOutlinks[i].second<<endl;
		}
		cout<<endl;
	}*/


}

bool gsjoin::sortbysec(const pair<int,double> &a,const pair<int,double> &b)
{
    return (a.second > b.second);
}
//ranks the vertices based on their pageranks and sorted
void gsjoin::rankingVertex()
{
	map<int,vector< pair<int,double> > >::iterator it;
	//map<int, vector< pair<int,int> > > sortedGraphs;
	for(it = idPageRank.begin(); it != idPageRank.end(); it++ )
	{
		int gid = it->first;
		vector< pair<int,double> > pagerankGraphs = it->second;
		sort(pagerankGraphs.begin(), pagerankGraphs.end(), sortbysec);
		/*cout << "After sorting vertices based on their pagerank" << endl;
		for(int i=0; i< pagerankGraphs.size(); i++){
			cout << pagerankGraphs[i].first << " " << pagerankGraphs[i].second << endl;
		}*/
		vector< pair<int,int> > vertexRanks;
		int rank = 1;
		for(int i = 0; i< pagerankGraphs.size(); i++)
		{
			vertexRanks.push_back(make_pair(pagerankGraphs[i].first, rank));
			rank++;
		}
		sortedGraphs[gid] = vertexRanks;
	}

	//printing sortedGraphs
	cout << "Graphs after ranking vertices based on their pagerank " << endl;
	for(map< int,vector< pair<int,int> > >::iterator it = sortedGraphs.begin();it != sortedGraphs.end();it++ ){
		cout << it->first << endl;
		vector< pair<int,int> > vertexRanks = it->second;
		for(int i=0; i< vertexRanks.size(); i++){
			cout << vertexRanks[i].first << " " <<vertexRanks[i].second << endl;
		}
		cout<<endl;
	}
	
}
void gsjoin::printGraphShingles(set<int>& graphShinglesSet)
{
	for(set<int>::iterator i = graphShinglesSet.begin(); i != graphShinglesSet.end(); i++)
		cout << *i << " ";
	cout << endl;
}
void gsjoin::printGraphShinglesBool(vector<int>& graphBool)
{
	for(vector<int>::iterator i = graphBool.begin(); i != graphBool.end(); i++)
		cout << *i << " ";
	cout << endl;
}
void gsjoin::printGraphSignature(vector<unsigned int>& graphSig)
{
	for(vector<unsigned int>::iterator i = graphSig.begin(); i != graphSig.end(); i++)
		cout << *i << " ";
	cout << endl;
}

void gsjoin::getShingles()
{
		map<int,vector< pair<int, int> > >::iterator it;
		size_t i;
		for(it = sortedGraphs.begin(); it != sortedGraphs.end(); it++)
		{
			set<int> shinglesSet;
			int gid = it->first;
			vector<pair<int, int> > verticesSequence = it->second;
			//cout << "DEBUG: getShingles1" << endl;
			for(i = 0; i < verticesSequence.size(); i++)
			{
				//cout << "DEBUG:Inside: getShingles3" << endl;
				shinglesSet.insert(verticesSequence[i].first);
			}
			graphShingles[gid] = shinglesSet;
			//cout << "DEBUG: getShingles2" << endl;
			//graphShingles[gid].insert(shinglesSet);
		}

}
void gsjoin::computeUnionShingles()
{
	map<int, set<int> >::iterator sit;
	for(sit = graphShingles.begin(); sit != graphShingles.end(); sit++)
	{
		set<int> shinglesSet = sit->second;
		unionShingles.insert(shinglesSet.begin(), shinglesSet.end());
		/*for(set<int>::iterator i = shinglesSet.begin(); i != shinglesSet.end(); i++)
		{
			unionShingles.insert(*i);
		}*/
	}
	cout << "All Shingles:";
	for(set<int>::iterator i = unionShingles.begin(); i != unionShingles.end(); i++)
		{
			cout << *i << ' ' << endl;
		}
}

void gsjoin::computeHashCoeffs()
 {
    // h = (a * x) + b
    // a and b should be randomly generated
    int i; 
    int MAX = numeric_limits<int>::max();
    for (i = 0; i < gsjoin::numHash; i++) {
     	hashCoeffs[i][0] = rand() % MAX +1; // a
        hashCoeffs[i][1] = rand() % MAX +1; // b
        cout << "DEBUG:hashCoeffs" << hashCoeffs[i][0] << ' ' << hashCoeffs[i][1] << endl;
    }
}

unsigned int gsjoin::computeHash(int index, int x, int numShingles) 
{
    return ((hashCoeffs[index][0] * x + hashCoeffs[index][1]) % numShingles);
}

int gsjoin::getMin(int a, int b)
{
	if(a < b)
		return a;
	return b;
}


  //returns Boolean vector of a graph for all shingles
vector<int> gsjoin::getBoolShingles(set<int>& vertexSet)
{
	set<int>::iterator it1,it2;
	//int i = 0;
	vector<int> graphShinglesBool;
	//cout << "DEBUG: getBoolShingles1" << endl;
	for(it1 = unionShingles.begin(); it1 != unionShingles.end(); it1++)
	{
		//cout << "DEBUG: getBoolShingles2" << endl;
		const bool is_in = vertexSet.find(*it1) != vertexSet.end();
		//vertexSet.count(*it1)
		if(is_in)
			graphShinglesBool.push_back(1);
		else
			graphShinglesBool.push_back(0);
		//cout << "DEBUG: getBoolShingles4" << endl;
	}
	return graphShinglesBool;
}
//returns vector of minHash Values for a graph
vector<unsigned int> gsjoin::getGraphMinHash(vector<int>& graphBool, int numShingles)
{
	set<int>::iterator it;
	int i;
	vector<unsigned int> graphSignatures;
	int MAX = numeric_limits<int>::max();
	for(i = 0;i < gsjoin::numHash; i++)
	{
		graphSignatures.push_back(MAX);
	}
	/*cout << "Signatures before minHash values" << endl;
	for(vector<unsigned int>::iterator it = graphSignatures.begin(); it != graphSignatures.end(); it++)
	{
		cout << *it << ' ';
	}*/
	//cout << endl;
	for(int r = 0; r < numShingles; r++)
	{
		if(graphBool[r] == 1)
		{
			for(i = 0; i < gsjoin::numHash; i++)
			{
				graphSignatures[i] = getMin(graphSignatures[i], computeHash(i,r,numShingles));
				//cout << graphSignatures[i] << ' ';
			}
			//cout << endl;
		}
	}
	return graphSignatures;
}

void gsjoin::performJoinSequenceSimilarity()
{
	map<int, set<int> >::iterator it1, it2;
	result resultFound;
	int numberOfPairs = 0;
	int numberofResults = 0;
	int numShingles = unionShingles.size();
	cout << "#Shingles: " << numShingles << endl;
	for(it1 = graphShingles.begin(); it1 != graphShingles.end(); it1++)
	{
		resultFound.xid = it1->first;
		set<int> graphShingles1 = graphShingles[resultFound.xid];
		vector<int> graphShinglesBool1 = getBoolShingles(graphShingles1);
		vector<unsigned int> graphSignatures1 = getGraphMinHash(graphShinglesBool1, numShingles);
		//printGraphShingles(graphShingles1);
		//printGraphShinglesBool(graphShinglesBool1);
		//printGraphSignature(graphSignatures1);
		for(it2 = graphShingles.begin(); it2 != graphShingles.end(); it2++)
		{
			resultFound.yid = it2->first;
			if(resultFound.xid == resultFound.yid || resultFound.xid > resultFound.yid)
				continue;
			set<int> graphShingles2 = graphShingles[resultFound.yid];
			vector<int> graphShinglesBool2 = getBoolShingles(graphShingles2);
			vector<unsigned int> graphSignatures2 = getGraphMinHash(graphShinglesBool2,numShingles);
			//printGraphShingles(graphShingles2);
			//printGraphShinglesBool(graphShinglesBool2);
			//printGraphSignature(graphSignatures2);
			double s = 0;
			for(int i = 0; i < graphSignatures1.size(); i++)
			{
				if(graphSignatures1[i] == graphSignatures2[i])
					s +=1;
			}
			resultFound.sim = s / graphSignatures1.size();
			 if(resultFound.sim >= gsjoin::threshold)
			{
				joinResult.push_back(resultFound);
				numberofResults++;
			}
			numberOfPairs++;
		}
	} 
}

//created by :jeevan
//to rank the vertices of by considering and pagerank value and adjacent vertices
void gsjoin::rankingVertexNew()
{
	//MODIFIED SIRS CODE TO RANK VERTICES BASED ON SORTING ORDER
		map<int,vector< pair<int,double> > >::iterator it;
		for(it = idPageRank.begin(); it != idPageRank.end(); it++ )
		{
			int gid = it->first;
			vector< pair<int,double> > pagerankGraphs = it->second;
			sort(pagerankGraphs.begin(), pagerankGraphs.end(), sortbysec);
			vector< pair<int,int> > vertexRanks;
			int rank = 1;
			for(int i = 0; i< pagerankGraphs.size(); i++)
			{
				vertexRanks.push_back(make_pair(pagerankGraphs[i].first, rank));
				rank++;
			}
			sortedGraphsOld[gid] = vertexRanks;
			getSortedSequence(gid);
		}	
		//printing sortedGraphs
		cout << "Graphs after ranking vertices based on their pagerank " << endl;
		for(map< int,vector< pair<int,int> > >::iterator it = sortedGraphs.begin();it != sortedGraphs.end();it++ ){
			cout << it->first << endl;
			vector< pair<int,int> > vertexRanks = it->second;
			for(int i=0; i< vertexRanks.size(); i++){
				cout << vertexRanks[i].first << " " <<vertexRanks[i].second << endl;
			}
			cout<<endl;
		}
}


//function to get the ranking sequnce based on the highest pagerank and adjacency
void  gsjoin::getSortedSequence(int gid)
{	
	
	vector< pair<int,int> > vertexRanksOldVector = sortedGraphsOld[gid];//to hold old vector
	map<int,int> vertexRanksOldMap; //map to hold the old vector map
	vector<int> visited; // to keep track of all the visited vertices (whose new page rank is found)
	if(vertexRanksOldVector.empty())
	{
		return;
	}
	int rank = 0; //to give the ranks
	vector< pair<int,int> >::iterator vertexRanksOldVecItr = vertexRanksOldVector.begin();
 	for(;vertexRanksOldVecItr != vertexRanksOldVector.end();vertexRanksOldVecItr++)
	{
		vertexRanksOldMap[vertexRanksOldVecItr->first] = vertexRanksOldVecItr->second;
	}
	vector< pair<int,int> > vertexRanks;//to hold new vector

	int CurrHighRankVertex = (vertexRanksOldVector.begin()->first);
	vertexRanks.push_back(make_pair(CurrHighRankVertex,++rank));//inserting into new vector
	 //deleting the from the old vector
	vertexRanksOldVector.erase(remove(vertexRanksOldVector.begin(),vertexRanksOldVector.end(),make_pair(CurrHighRankVertex,vertexRanksOldMap[CurrHighRankVertex])),vertexRanksOldVector.end());
	visited.push_back(CurrHighRankVertex);
	vector<int> adjList;
	while(!vertexRanksOldVector.empty())
	{
		adjList = getAdjList(gid,CurrHighRankVertex); //get adjacent list of current highest page ranked vertex
		if(!adjList.empty()) //if adjcennt list is not empty
		{	
			while(!adjList.empty() && find(visited.begin(),visited.end(),adjList.front()) != visited.end())//if the vertex is already visited
			{
				adjList.erase(adjList.begin()); //delete that vertex from  adjacency list
			}
			if(!adjList.empty())
			{
				//get the highest element from the vector and insert into vertexrank
				CurrHighRankVertex = adjList.front();
				adjList.erase(adjList.begin());//delete the vertex  from adjacency list
				vertexRanks.push_back(make_pair(CurrHighRankVertex,++rank)); //insert to new vector
				//deleting the vertex from the old vector
				vertexRanksOldVector.erase(remove(vertexRanksOldVector.begin(),vertexRanksOldVector.end(),make_pair(CurrHighRankVertex,vertexRanksOldMap[CurrHighRankVertex])),vertexRanksOldVector.end());
				visited.push_back(CurrHighRankVertex);
				continue;
			}
		}
		CurrHighRankVertex = (vertexRanksOldVector.begin())->first;
		vertexRanks.push_back(make_pair(CurrHighRankVertex,++rank)); //insert to new vector
		//deleting the from the old vector
		vertexRanksOldVector.erase(remove(vertexRanksOldVector.begin(),vertexRanksOldVector.end(),make_pair(CurrHighRankVertex,vertexRanksOldMap[CurrHighRankVertex])),vertexRanksOldVector.end());
		visited.push_back(CurrHighRankVertex);
	}
	//inserting new vector of vertexranks onto sorted graphs
	sortedGraphs[gid] = vertexRanks;
}

//function to get the adjacent list from given graph id and vertex id
vector<int> gsjoin::getAdjList(int gid,int vertex){
	vector< pair<int,int> > EdgeList;
	//get the graph representation(EdgeList) for the given graph id from idGraph Map
	map<int,vector< pair<int,int> > >::iterator idGpIt;
	for(idGpIt = idGraph.begin();idGpIt != idGraph.end();idGpIt++)
	{
		if(idGpIt->first == gid)
		{
			EdgeList = idGpIt->second;
			break;
		}
	}
	vector<int> adjVertices;
	vector< pair<int,int> >::iterator EdgeListItr;
	//get the adjacent vertices list for the given vertex
	for(EdgeListItr = EdgeList.begin();EdgeListItr != EdgeList.end() ; EdgeListItr++)
	{
		if(EdgeListItr->first == vertex)
		{
			adjVertices.push_back(EdgeListItr->second);
		}
	}
	//the adjacent vertices list may contain the duplicate entries 
	// remove the duplicate entries using unique function
	vector<int>::iterator 	 AdjVerticesItr;
	sort(adjVertices.begin(),adjVertices.end());
	adjVertices.erase( unique( adjVertices.begin(), adjVertices.end() ), adjVertices.end() );
	//create a vector of pairs (vertex,pagerank)
		//1.get the pagerank vector for given graph
		vector< pair<int,double> >pageRanks = idPageRank[gid];
		//2.create the vector of pairs
		vector< pair<int,double> >::iterator pageRanksItr; //iterator to traverse pageRanks
		//3.insert the pageranks into map for better accessiblity
			map <int,double>  pageRanksMap;
			for(pageRanksItr = pageRanks.begin();pageRanksItr != pageRanks.end();pageRanksItr++)
			{
				pageRanksMap[pageRanksItr->first] = pageRanksItr->second;
			}
			vector< pair<int,double> > outputPageRanks;
			for(AdjVerticesItr = adjVertices.begin();AdjVerticesItr != adjVertices.end();AdjVerticesItr++)
			{
				outputPageRanks.push_back(make_pair(*AdjVerticesItr, pageRanksMap[*AdjVerticesItr]));
			}
		//4.sort the output pageranks vector based on the pageranks
		sort(outputPageRanks.begin(),outputPageRanks.end(),sortbysec);
		//5.extract the sorted elements into new output vector
		vector<int> outputVec;
		for(pageRanksItr = outputPageRanks.begin();pageRanksItr != outputPageRanks.end();pageRanksItr++)
		{
			outputVec.push_back(pageRanksItr->first);
		}			
	return outputVec;
}
//--------- end of sequence similarity ----------------
//author: jeev@n



// -------------  start of constructing inverted index and comparing the graphs using inverted index -----------------
// author : jeev@n

//to construct inverted index for given graphs
void gsjoin::constructInvertedIndex()
{
	set<int>::iterator unionShinglesItr;
	map<int , set<int> >::iterator graphShinglesItr;
	set<int> CurrentGraphSingles; 
	set<int>::iterator CurrentGraphSinglesItr;
	int currentGID;
	set<int> lst; // to store the list of GID's for given shingle
	set<int>::iterator listItr; // to iterate list of GID's
	set<int>::iterator listitr1; 
	
	//constructing inverted index from union shingles and graph shingles
	cout << endl << "------ Inverted Index -------" ; 
	//for eve\ry shingle in the set of union shingles
	for(unionShinglesItr = unionShingles.begin(); unionShinglesItr != unionShingles.end();unionShinglesItr++)
	{
		//construsting list for current shingle
		cout << endl;
		cout << *(unionShinglesItr) << " : " ;
		//check how many graphs contains this shingle in it
		for(graphShinglesItr = graphShingles.begin();graphShinglesItr != graphShingles.end();graphShinglesItr++)
		{
			 // cout << "error1" << endl;
			currentGID = graphShinglesItr->first;
			CurrentGraphSingles = graphShinglesItr->second;
			for(CurrentGraphSinglesItr = CurrentGraphSingles.begin();CurrentGraphSinglesItr != CurrentGraphSingles.end();CurrentGraphSinglesItr++)
			{
				//if a graph contains this shingle 
				if((*unionShinglesItr) == (*CurrentGraphSinglesItr))
				{
					//add the graph id in inverted index
					invertedIndex[*unionShinglesItr].insert(currentGID);
					cout << currentGID  << " -> " ; 
				}
			}
		}
		//constructing related vertices using aboeve list
		lst = invertedIndex[*unionShinglesItr];
		for(listItr = lst.begin(); listItr != prev(lst.end());listItr++)
		{
			listitr1 = listItr;
			listitr1++;
			for(;listitr1 != lst.end();listitr1++)
			{
				relatedGaphs[*listItr].insert(*listitr1);
			}
		}
	}
	//printing related graphs
	cout << endl << "----  printing related GRAPHS ---------";
  	map<int, set<int> >::iterator relatedGaphsItr;
	for(relatedGaphsItr = relatedGaphs.begin();relatedGaphsItr != relatedGaphs.end();relatedGaphsItr++)
	{
		cout << endl;
		cout << relatedGaphsItr->first << " : " ;
		lst = relatedGaphsItr->second;
		for(listItr = lst.begin();listItr != lst.end();listItr++)
		{
			cout << *listItr << " -> " ;
		}
	}
	cout << endl << endl;
}

//to perform join sequence smilarilty by using inverted index
void gsjoin::performJoinSequenceSimilarityNew()
{

	result resultFound; //to store the results
	set<int> graphSet; 
	set<int>::iterator graphSetItr;
	int numberOfPairs = 0;
	int numberofResults = 0;
	int numShingles = unionShingles.size();
	cout << "#Shingles: " << numShingles << endl;
	map<int,set<int> >::iterator relatedGaphsItr;
	//traversing related graphs 
	for(relatedGaphsItr = relatedGaphs.begin();relatedGaphsItr != relatedGaphs.end();relatedGaphsItr++)
	{
		resultFound.xid = relatedGaphsItr->first;
		graphSet = relatedGaphsItr->second;
		//get the signature for the first graph
		set<int> graphShingles1 = graphShingles[resultFound.xid];
		vector<int> graphShinglesBool1 = getBoolShingles(graphShingles1);
		vector<unsigned int> graphSignatures1 = getGraphMinHash(graphShinglesBool1, numShingles);
		//traversing graphID set for every GID in related graphs
		for(graphSetItr = graphSet.begin();graphSetItr != graphSet.end(); graphSetItr++)
		{
			resultFound.yid = *graphSetItr;
			//if it is already computed or not necessary then continue
			if(resultFound.xid == resultFound.yid || resultFound.xid > resultFound.yid)
				continue;
			//compute the graph signaute for second graph
			set<int> graphShingles2 = graphShingles[resultFound.yid];
			vector<int> graphShinglesBool2 = getBoolShingles(graphShingles2);
			vector<unsigned int> graphSignatures2 = getGraphMinHash(graphShinglesBool2,numShingles);
			//compare the signatures and compute the similarity among them
			double s = 0;
			for(int i = 0; i < graphSignatures1.size(); i++)
			{
				if(graphSignatures1[i] == graphSignatures2[i])
					s +=1;
			}
			resultFound.sim = s / graphSignatures1.size();
			//if the similarity found is greater than threshold then save the result
			 if(resultFound.sim >= gsjoin::threshold)
			{
				joinResult.push_back(resultFound);
				numberofResults++;
			}
			numberOfPairs++;
		}
	}
}
//------- end of inverted index and comparing graphs using it-------------------------------------------------------
