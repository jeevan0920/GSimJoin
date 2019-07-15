#include <sys/time.h>
#include <ctime>
#include <fstream>

//#include "gsjoin.h"
#include "gsjoin.cpp"
#include "io.cpp"

int main(int argc, char* argv[])
{
	gsjoin gsjoinInstance;
	if (argc < 4) {
		cerr << "Usage: gsjoin [sim_func] [threshold] [input_file]" << endl;
		return 1;
	}

	else {
		if (gsjoinInstance.ReadParameters(argv[1], argv[2])) {
			cerr << "Similarity Function Missing!" << endl;
			return 1;
		}
	}

	cerr << "Document: " << argv[3] << endl;

	
	gsjoinInstance.readDataNew(argv[3]);
	cout << "successs" << endl;
	gsjoinInstance.getPageRanks();

	
	clock_t c_start, c_end;
	c_start = clock();	
	timeval timeStart, timeEnd;
	gettimeofday(&timeStart, NULL);
	cerr << "=== BEGIN JOIN (TIMER STARTED) ===" << endl;
	//gsjoinInstance.rankingVertex();
	gsjoinInstance.rankingVertexNew();
	gsjoinInstance.getShingles();
	gsjoinInstance.computeUnionShingles();
	gsjoinInstance.computeHashCoeffs();
	//gsjoinInstance.VertexRanking();
	gsjoinInstance.constructInvertedIndex();
	gsjoinInstance.performJoinSequenceSimilarityNew();
	gsjoinInstance.WriteResults();
	

	cerr << "=== END JOIN (TIMER STOPPED) ===" << endl;
	c_end = clock();
	gettimeofday(&timeEnd, NULL);
	//clock_gettime(CLOCK_REALTIME, &timeEnd);
	cerr << "Total Running Time: " << setiosflags(ios::fixed) << setprecision(6) << double(timeEnd.tv_sec - timeStart.tv_sec) + double(timeEnd.tv_usec - timeStart.tv_usec) / 1e6 << " sec" << endl;
	//cerr << "Total Running Time: " << setiosflags(ios::fixed) << setprecision(6) << double(timeEnd.tv_sec - timeStart.tv_sec) + double(timeEnd.tv_nsec - timeStart.tv_nsec) / 1e9 << endl;
	cerr << "Total Running Time: " << setiosflags(ios::fixed) << setprecision(3) << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms" << endl;
	return 0;
}
