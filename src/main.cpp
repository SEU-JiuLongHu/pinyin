#include <iostream>
#include <fstream>
#include "HMM.h"
using namespace std;

int main() {

	HMM hmm;
	//system("dir");
	ifstream fi("in.txt");
	string tmp;
	fi >> tmp;
	if (hmm.loadPinyinCharMap("F:\\Github\\pinyin\\py2word.txt"))
	{
		if (hmm.loadCorpus("F:\\Github\\pinyin\\corpus.txt"))
		{
			hmm.query(vector<string>{"mei","you","an","zhuang","qiao"},10);
			cout << "done" << endl;
		}
	}

	//cout << numeric_limits<double>().lowest() << endl;
	system("pause");
	return 0;
}