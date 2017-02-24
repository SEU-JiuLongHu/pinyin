#include <iostream>
#include "../include/common.hpp"
#include "../include/pyt.hpp"

using namespace std;
using namespace pinyin;

int main() {
	PinYinTrie pinYinTrie;
	pinYinTrie.Init();
	pinYinTrie.Build();
	std::vector<std::string> res = pinYinTrie.SplitPinYin("tiananmen");
	for (std::string t : res) {
		std::cout << t << std::endl;
	}
	system("pause");
	return 0;
}