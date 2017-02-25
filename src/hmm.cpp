#include <iostream>
#include "../include/hmm.hpp"
#include "../include/common.hpp"

using namespace std;

namespace pinyin {

HMM::HMM() {
	std::cout << "" << std::endl;
}

std::map<std::string, double> HMM::PY2Chinese(std::string pys) {
	vector<std::string> splitPYs = this->pyTrie.SplitPinYin(pys);
	std::map<std::string, double> chineseRes;
	for (std::string py : splitPYs) {
		std::map<std::string, double> _r = Viterbi(py);
		chineseRes.insert(_r.begin(), _r.end());
	}
	return chineseRes;
}

std::map<std::string, double> HMM::Viterbi(std::string py) {
	std::vector<std::string> elems = split(py, '+');
	std::map<std::string, double> V = this->hmmTable.QueryStarting(elems[0]);
	for (size_t i = 1; i < elems.size(); i++) {
		std::string _py = elems[i];
		std::map<std::string, double> probM;
		for (auto const& pp: V) {
			std::string character = pp.first;
			double prob = pp.second;
			std::map<std::string, double> rr = this->hmmTable.QueryTransfer(_py, character);
			if (0 == rr.size()) {
				continue;
			}
			for (auto const& _p : rr) {
				probM[character + _p.first] = _p.second + prob;
			}
		}
		if (0 == probM.size()) {
			return V;
		} else {
			V = probM;
		}
	}
	return V;
}
}