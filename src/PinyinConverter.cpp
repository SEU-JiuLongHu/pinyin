#include "PinyinConverter.h"
PinyinConverter::PinyinConverter()
{
}

PinyinConverter::~PinyinConverter()
{
}

string PinyinConverter::UnicodeToUtf8(wstring unicode)
{
	static wstring_convert<codecvt_utf8<wchar_t>> converter;
	string utf8Str = converter.to_bytes(unicode);
	return utf8Str;
}

wstring PinyinConverter::Utf8ToUnicode(string utf8)
{
	static wstring_convert<codecvt_utf8<wchar_t>> converter;
	wstring unicodeStr = converter.from_bytes(utf8);
	return unicodeStr;
}
