#include <iostream>
#include <vector>
#include <string.h>
using namespace std;

vector<char*> raw_txt;
int stringNum;

void ReadRawText(char* filename)
{
	char* inputString = new char [1000000];
	stringNum = 0;
	freopen(filename, "r", stdin);
	while (gets(inputString)) {
		char* newString = new char [strlen(inputString) + 1];
		strcpy(newString, inputString);
		raw_txt.push_back(newString);
		++stringNum;
	}
	delete inputString;
}

void ReadResult(char* filename)
{
	int x_id, y_id;
	double sim;
	int resultNum = 1;
	freopen(filename, "r", stdin);
	while (scanf("%d %d %lf", &x_id, &y_id, &sim) == 3) {
		printf("Result Pair %d:\n", resultNum);
		printf("%s\n%s\n", raw_txt[x_id - 1], raw_txt[y_id - 1]);
		printf("Similarity: %.3f\n\n", sim);
		++resultNum;
	}
}

void FreeRawText()
{
	for (int i = 0; i < raw_txt.size(); i++) 
		delete [] raw_txt[i];
}

int main(int argc, char* argv[])
{
	if (argc < 3) {
		cerr << "Usage: showresult [raw_txt_file] [result_file]" << endl;
		return 1;
	}

	ReadRawText(argv[1]);
	ReadResult(argv[2]);
	FreeRawText();

	return 0;
}
