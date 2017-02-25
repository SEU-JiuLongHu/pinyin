#include <iostream>
#include <fstream>
#include "../include/common.hpp"
#include "../include/hmm.hpp"

using namespace std;
using namespace pinyin;

int main() {
	/*
	PinYinTrie pinYinTrie;
	std::vector<std::string> res = pinYinTrie.SplitPinYin("tiananmen");
	for (std::string t : res) {
		std::cout << t << std::endl;
	}
	*/
	
	//ifstream cin("in.txt");
	ofstream cout("out.txt");
	//HMMTable hmmtable;
	
	//std::map<std::string, double> res = hmmtable.QueryStarting("xie");	
	
	//std::wstring _wstr ("");
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
	std::wstring wpre(L"°Ö");
	
	std::string pre = converter.to_bytes(wpre);
	//std::map<std::string, double> res = hmmtable.QueryTransfer("ma", pre);

	/*
	std::map<std::string, double> r;
	r["a"] = 1;
	r["b"] = 2;
	std::map<std::string, double> c;
	c["c"] = 4;
	c["d"] = 5;
	r = c;
	//r.insert(c.begin(), c.end());
	for (auto const& pp : r) {
		cout << pp.first << " " << pp.second << endl;
	}
	*/

	HMM hmm;
	std::vector<std::pair<std::string, double>> res = hmm.PY2Chinese("wohenkuaile");
	for (auto const& pp: res) {
		//std::string narrow = converter.to_bytes(wide_utf16_source_string);
		//std::wstring wide = converter.from_bytes(res[i][j]);
		//std::string t = PinyinConverter::UnicodeToUtf8(wide);
		cout << pp.first << " " << pp.second << endl;
	}
	cout.close();
	system("pause");
	return 0;
}