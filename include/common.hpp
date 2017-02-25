#include <codecvt>
#include <locale>
#include <string>
#include <vector>

using std::string;
using std::wstring;
using std::wstring_convert;
using std::codecvt_utf8;
using namespace std;

namespace pinyin {

class Config {
public:
	static inline Config& GetInstance() {
		static Config c;
		return c;
	}

public:
	std::string dataset;
	char* hhmdbname;

private:
	Config();
	~Config() {}
	Config(const Config &other);
	Config &operator=(const Config &other);
};

class PinyinConverter {
public:
	PinyinConverter();
	~PinyinConverter();
	static string UnicodeToUtf8(wstring unicode);
	static wstring Utf8ToUnicode(string utf8);
};

void split(const std::string &, char , std::vector<std::string>& );
std::vector<std::string> split(const std::string &, char );
}