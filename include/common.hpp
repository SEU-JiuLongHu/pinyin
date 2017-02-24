#include <string>

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

private:
	Config();
	~Config() {}
	Config(const Config &other);
	Config &operator=(const Config &other);
};
}