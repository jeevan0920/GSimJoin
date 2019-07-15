//#include "gsjoin.h"
#include <algorithm>

void gsjoin::CalcPrefix()
{
	int i, j, len;
	int widowTokenNum = 0;
	int indexPrefixTokenNum = 0;
	tokenFrequency.resize(tokenNum);
	for (i = 0; i < tokenNum; i++) tokenFrequency[i] = 0; // initialize token frequency

	for (i = 0; i < recordNum; i++) {
		len = joinRecord[i].len;

		/* calculate indexing and probing prefix lengths */
		if (jaccard) {
			//joinRecord[i].indexPrefixLength = (int)(len - len * 2.0 * threshold / (1.0 + threshold) + EPS) + 1;
			joinRecord[i].indexPrefixLength = (int)(len - len * threshold + EPS) + 1;
			joinRecord[i].probePrefixLength = (int)(len - len * threshold + EPS) + 1;
		}
		else if (cosine) {
			joinRecord[i].indexPrefixLength = (int)(len - len * threshold + EPS) + 1;
			joinRecord[i].probePrefixLength = (int)(len - len * threshold * threshold + EPS) + 1;
		}

		if (joinRecord[i].indexPrefixLength > len) joinRecord[i].indexPrefixLength = len;
		if (joinRecord[i].probePrefixLength > len) joinRecord[i].probePrefixLength = len;

		/* count token frequency */
		for (j = 0; j < joinRecord[i].probePrefixLength; j++)
			++tokenFrequency[joinRecord[i].tokenArray[j]];
	}

	/* count widow tokens and non-widow tokens */
	for (i = 0; i < tokenNum; i++) {
		if (tokenFrequency[i] == maxWidowTokenFrequency) ++widowTokenNum;
		if (tokenFrequency[i]) ++indexPrefixTokenNum;
	}
	cerr << "# Distinct Tokens: " << tokenNum << endl;
	cerr << "# Distinct Indexable Tokens: " << indexPrefixTokenNum << endl;//max #of tokens can be indexed
	cerr << "# Distinct Widow Tokens: " << widowTokenNum << endl;//Tokens with frequency one in prefixLength
}

void gsjoin::BuildIndex()
{
	int i, j, len, token;
	int indexEntryNum = 0;
	index_entry indexEntry;

	/* indexing phase */
	invertedList.resize(tokenNum);
	for (i = 0; i < tokenNum; i++) {
		invertedList[i].startingPosition = -1;
	}
	for (i = 0; i < recordNum; i++) {
		len = joinRecord[i].len; 
		for (j = 0; j < joinRecord[i].indexPrefixLength; j++) {
			token = joinRecord[i].tokenArray[j]; // fetch token

			if (tokenFrequency[token] <= maxWidowTokenFrequency) continue; // widow token

			/* index record id and its position */
			indexEntry.rec = i, indexEntry.pos = j;
			invertedList[token].startingPosition = 0;
			invertedList[token].list.push_back(indexEntry);

			++indexEntryNum;
		}
	}

	/* output statistics */
	cerr << "# Inverted Index Entries: " << indexEntryNum << endl;
}

