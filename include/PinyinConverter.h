#pragma once
#include <codecvt>
#include <locale>
#include <string>
using std::string;
using std::wstring;
using std::wstring_convert;
using std::codecvt_utf8;
class PinyinConverter
{
public:
	PinyinConverter();
	~PinyinConverter();
	static string UnicodeToUtf8(wstring unicode);
	static wstring Utf8ToUnicode(string utf8);
};

