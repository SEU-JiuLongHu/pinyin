#include <iostream>
#include "../include/common.hpp"
#include "../include/pyt.hpp"
#include "../include/hmmdb.hpp"
#include <fstream>
using namespace std;
using namespace pinyin;

int main() {
	/*
	PinYinTrie pinYinTrie;
	pinYinTrie.Init();
	pinYinTrie.Build();
	std::vector<std::string> res = pinYinTrie.SplitPinYin("tiananmen");
	for (std::string t : res) {
		std::cout << t << std::endl;
	}
	*/
	ifstream cin("in.txt");
	ofstream cout("out.txt");
	HMMTable hmmtable;
	
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
	std::wstring wpre(L"Íõ");
	std::string pre = converter.to_bytes(wpre);
	
	//std::vector<std::vector<std::string>> res = hmmtable.QueryStarting("xie");
	//std::wstring _wstr ("");
	std::vector<std::vector<std::string>> res = hmmtable.QueryTransfer("ba", pre);
	for (size_t i = 0; i < res.size(); i++) {
		for (size_t j = 0; j < res[i].size(); j++) {
			//std::string narrow = converter.to_bytes(wide_utf16_source_string);
			//std::wstring wide = converter.from_bytes(res[i][j]);
			//std::string t = PinyinConverter::UnicodeToUtf8(wide);
			cout << res[i][j].c_str() << " ";
		}
		cout << std::endl;
	}
	system("pause");
	return 0;
}