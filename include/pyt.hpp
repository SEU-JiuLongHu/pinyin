#include <string>
#include <vector>

using namespace std;

#define CHAR_RANGE 26

namespace pinyin {

struct TrieNode {
	bool isPY;
	TrieNode* parent;
	TrieNode* childs[CHAR_RANGE];
	std::string py;
	TrieNode() : parent(NULL), isPY(false) {
		for (size_t i = 0; i < CHAR_RANGE; i++) {
			childs[i] = NULL;
		}
	}
};

class PinYinTrie {
public:
	PinYinTrie();
	~PinYinTrie();
	void Init();
	void Build();
	std::vector<std::string> SplitPinYin(std::string);

public:
	std::string dataset;
	TrieNode* root;

private:
	std::vector<string> ReadPinYins();
	void InsertDFS(TrieNode* , std::string, int );
	void SearchDFS(TrieNode*, int, std::string, vector<std::string> & );
	void BuildTrie(std::vector<std::string> & );
};
}