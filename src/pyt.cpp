#include <iostream>
#include <fstream>
#include <sstream>
#include "../include/pyt.hpp"
#include "../include/common.hpp"

namespace pinyin {

PinYinTrie::PinYinTrie() {
	this->Init();
	this->Build();
}
PinYinTrie::~PinYinTrie() {}

void PinYinTrie::Init() {
	Config& config = Config::GetInstance();
	this->dataset = config.dataset;
}

std::vector<std::string> PinYinTrie::ReadPinYins() {
	std::vector<std::string> pinyins;
	std::ifstream infile(this->dataset);
	std::string line;
	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		std::string p;
		iss >> p;
		//std::cout << p << std::endl;
		pinyins.push_back(p);
	}
	return pinyins;
}

void PinYinTrie::InsertDFS(TrieNode* cur, std::string word, int idx) {
	if (idx == word.size()) {
		cur->isPY = true;
		cur->py = word;
		return;
	}
	int index = word[idx] - 'a';
	if (index < 0 || index >= 26) {
		std::cout << word[idx] << std::endl;
		std::cout << "Unknow PinYin Char" << std::endl;
	}
	if (NULL == cur->childs[index]) {
		cur->childs[index] = new TrieNode();
	}
	InsertDFS(cur->childs[index], word, idx + 1);
}

void PinYinTrie::BuildTrie(std::vector<std::string>& pinyins) {
	this->root = new TrieNode();
	for (size_t i = 0; i < pinyins.size(); i++) {
		InsertDFS(root, pinyins[i], 0);
	}
}

void PinYinTrie::Build() {
	std::vector<std::string> pinyins = ReadPinYins();
	BuildTrie(pinyins);
}

void PinYinTrie::SearchDFS(TrieNode* cur, int idx, std::string _pys, vector<std::string> &pres) {
	if (idx == _pys.size()) {
		if (cur->isPY) {
			pres.push_back(cur->py);
		}
		return;
	}
	if (cur->isPY) {
		pres.push_back(cur->py);
	}
	int index = _pys[idx] - 'a';
	if (NULL == cur->childs[index]) {
		return;
	}
	return SearchDFS(cur->childs[index], idx + 1, _pys, pres);
}

std::vector<std::string> PinYinTrie::SplitPinYin(std::string pys) {
	std::vector<std::string> res;
	std::string _pys = pys;
	std::vector<std::string> pres;
	SearchDFS(this->root, 0, _pys, pres);
	for (std::string t : pres) {
		int _l = t.size();
		if (_l == _pys.size()) {
			res.push_back(t);
			continue;
		}
		std::string leftpys = _pys.substr(_l, _pys.size());
		std::vector<std::string> leftpres = SplitPinYin(leftpys);
		
		for (std::string l : leftpres) {
			std::string _r = t + "+" + l;
			res.push_back(_r);
		}
	}
	return res;
}
}