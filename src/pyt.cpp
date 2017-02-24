#include <iostream>
#include <fstream>
#include <sstream>
#include "../include/pyt.hpp"
#include "../include/common.hpp"

namespace pinyin {

PinYinTrie::PinYinTrie() {}
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

void PinYinTrie::BuildTrie() {

}
}