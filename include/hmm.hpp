#include <string>
#include <vector>
#include "hmmdb.hpp"
#include "pyt.hpp"

using namespace std;

namespace pinyin {

class HMM {
public:
	HMM();
	std::map<std::string, double> PY2Chinese(std::string);

private:
	std::map<std::string, double> Viterbi(std::string);

private:
	HMMTable hmmTable;
	PinYinTrie pyTrie;
};
}