#include <iostream>
#include "../include/common.hpp"
#include "../include/pyt.hpp"

using namespace std;
using namespace pinyin;

int main() {
	PinYinTrie pinYinTrie;
	pinYinTrie.Init();
	pinYinTrie.ReadPinYins();
	system("pause");
	return 0;
}