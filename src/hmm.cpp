#include <iostream>
#include <algorithm>
#include "../include/hmm.hpp"
#include "../include/common.hpp"

using namespace std;

namespace pinyin {

HMM::HMM() {
	//std::cout << "" << std::endl;
}

std::vector<std::pair<std::string, double>> HMM::PY2Chinese(std::string pys) {
	vector<std::string> splitPYs = this->pyTrie.SplitPinYin(pys);
	std::map<std::string, double> chineseRes;
	for (std::string py : splitPYs) {
		std::map<std::string, double> _r = Viterbi(py);
		chineseRes.insert(_r.begin(), _r.end());
	}
	std::vector<std::pair<std::string, double>> cR;
	for (auto iter = chineseRes.begin(); iter != chineseRes.end(); iter++) {
		cR.push_back(*iter);
	}
	auto cmp = [](std::pair<std::string,double> const& a, std::pair<std::string, double> const& b) {
		return a.second != b.second ? a.second > b.second : a.first < b.first;
	};
	std::sort(cR.begin(), cR.end(), cmp);
	return cR;
}

std::map<std::string, double> HMM::Viterbi(std::string py) {
	std::vector<std::string> elems = split(py, '+');
	std::map<std::string, double> V = this->hmmTable.QueryStarting(elems[0]);
	for (size_t i = 1; i < elems.size(); i++) {
		std::string _py = elems[i];
		std::map<std::string, double> probM;
		for (auto const& pp: V) {
			std::string cc = pp.first;
			std::string character(cc.substr(cc.size() - 3, 3));
			double prob = pp.second;
			std::map<std::string, double> rr = this->hmmTable.QueryTransfer(_py, character);
			if (0 == rr.size()) {
				continue;
			}
			for (auto const& _p : rr) {
				probM[pp.first + _p.first] = _p.second + prob;
			}
		}
		if (0 == probM.size()) {
			// if no transfer
			// restart from current pinyin
			std::string _py = elems[i];
			for (size_t j = i + 1; j < elems.size(); j++) {
				_py += "+";
				_py += elems[j];
			}
			std::map<std::string, double> rV = Viterbi(_py);
			std::map<std::string, double> finalV;
			for (auto const& _pp : V) {
				for (auto const& rpp : rV) {
					finalV[_pp.first + rpp.first] = _pp.second + rpp.second;
				}
			}
			return finalV;
		} else {
			V = probM;
		}
	}
	return V;
}
}