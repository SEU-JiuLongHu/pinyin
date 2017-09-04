#include <iostream>
#include <fstream>
#include "HMM.h"
using namespace std;
struct pRes
{
	wchar_t *str;
	double score;
};
extern "C"
{
	HMM hmm;
	bool load_corpus(char *corpusname)
	{
		string filename(corpusname);
		return hmm.loadCorpus(filename);
	}
	pRes query(char **pinyins, int length)
	{
		vector<string> strs(length);
		for (size_t i = 0; i < length; i++)
			strs[i] = string(pinyins[i]);
		vector<Res> res = hmm.query(strs, 1);
		Res r = res[0];
		pRes pres;
		pres.score = r.prob;
		//int len = (r.str.c_str());
		pres.str = new wchar_t[r.str.length()+1];
		for (size_t i = 0; i < r.str.length(); i++)
			pres.str[i] = r.str[i];
		pres.str[r.str.length()] = '\0';
		return pres;
	}
}
int main() {

	//HMM hmm;
	//system("dir");
		//if (hmm.loadCorpus("F:\\Github\\pinyin\\corpus.txt"))
		//{
		//	hmm.query(vector<string>{"jiao","cheng"},10);
		//	cout << "done" << endl;
		//}

	//cout << numeric_limits<double>().lowest() << endl;
	load_corpus("F:\\Github\\pinyin\\corpus.txt");
	char *str[] = { "yan","jiu","sheng",NULL };
	pRes res = query(str, 3);
	int test = 0;
	test++;
	wstring tmp= wstring(res.str);
	
	int n = tmp.length();
	return 0;
}