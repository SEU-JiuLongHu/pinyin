#include <string>
#include <vector>
#include <map>
#include "../3rdparty/sqlite3.h"

using namespace std;

namespace pinyin {

class HMMTable {
public:
	HMMTable();
	HMMTable& operator=(const HMMTable);
	~HMMTable();
	
	std::map<std::string, double> QueryStarting(std::string);
	std::map<std::string, double> QueryTransfer(std::string, std::string);
	std::map<std::string, double> QueryObservation();
	
private:
	void ConnectDB();
	void DisConnectDB();

public:
	char* hmmdab;

private:
	sqlite3* pDB = NULL;
};
}