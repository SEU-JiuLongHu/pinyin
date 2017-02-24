#include <string>
#include <vector>

using namespace std;

namespace pinyin {

class PinYinTrie {
public:
	PinYinTrie();
	~PinYinTrie();

	void Init();
	std::vector<string> ReadPinYins();
	void Write();
	void BuildTrie();

public:
	std::string dataset;
};
}