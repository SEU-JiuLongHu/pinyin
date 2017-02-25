#include <string>
#include <sstream>
#include <vector>
#include "../include/common.hpp"

namespace pinyin {

Config::Config() {
	dataset = "../py2word.txt";
	hhmdbname = "../hmm.sqlite";
}

PinyinConverter::PinyinConverter() {}

PinyinConverter::~PinyinConverter() {}

string PinyinConverter::UnicodeToUtf8(wstring unicode) {
	static wstring_convert<codecvt_utf8<wchar_t>> converter;
	string utf8Str = converter.to_bytes(unicode);
	return utf8Str;
}

wstring PinyinConverter::Utf8ToUnicode(string utf8) {
	static wstring_convert<codecvt_utf8<wchar_t>> converter;
	wstring unicodeStr = converter.from_bytes(utf8);
	return unicodeStr;
}

void split(const std::string &s, char delim, std::vector<std::string>& result) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		result.push_back(item);
	}
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}
}