#include <iostream>
#include <string>
#include <vector>
#include "../include/hmmdb.hpp"
#include "../include/common.hpp"

using namespace std;

namespace pinyin {

HMMTable::HMMTable() {
	this->ConnectDB();
}

HMMTable& HMMTable::operator=(const HMMTable h) {
	this->pDB = h.pDB;
	return *this;
}

HMMTable::~HMMTable() {
	this->DisConnectDB();
}

void HMMTable::ConnectDB() {
	Config& config = Config::GetInstance();
	this->hmmdab = config.hhmdbname;
	int nRes = sqlite3_open(this->hmmdab, &this->pDB);
	if (nRes != SQLITE_OK) {
		std::cout << "Open database error." << std::endl;
		exit(0);
	}
	std::cout << "Open database success." << std::endl;
}

void HMMTable::DisConnectDB() {
	sqlite3_close(this->pDB);
	std::cout << "Close database success." << std::endl;
}

/*vector<string> HMMTable::SqlCallback(void *data, int argc, char **argv, char **azColName){
	vector<> res;
	for (int i = 0; i < argc, i++) {

	}
}*/

std::map<std::string, double> HMMTable::QueryStarting(std::string py) {
	const char *zTail;
	sqlite3_stmt* stmt = NULL;
	//std::vector<std::vector<std::string>> res;
	std::map<std::string, double> res;
	//std::string sql = "SELECT e.character, (e.probability + s.probability) FROM EMISSION e, STARTING s WHERE e.character==s.character AND e.pinyin='zui' ORDER BY (e.probability + s.probability) DESC LIMIT 10;";
	std::string sql = "SELECT e.character, (e.probability + s.probability) FROM EMISSION e, STARTING s WHERE e.character==s.character AND e.pinyin='" + py + "' ORDER BY (e.probability + s.probability) DESC LIMIT 10;";
	//std::cout << sql << std::endl;
	const char* _sql = sql.c_str();
	if (sqlite3_prepare_v2(this->pDB, _sql, -1, &stmt, &zTail) == SQLITE_OK) {
		while (sqlite3_step(stmt) == SQLITE_ROW) {
			//std::vector<std::string> r;
			std::string _py(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0)));
			double pro = sqlite3_column_double(stmt, 1);
			//r.push_back(_py);
			//r.push_back(_prob);
			//res.push_back(r);
			res[_py] = pro;
		}
	}
	return res;
}

std::map<std::string, double> HMMTable::QueryTransfer(std::string py, std::string chinese) {
	// chinese: former chinese
	const char *zTail;
	sqlite3_stmt* stmt = NULL;
	//std::vector<std::vector<std::string>> res;
	std::map<std::string, double> res;
	std::string sql = "SELECT t.behind, (e.probability + t.probability) FROM EMISSION e, TRANSITION t WHERE e.character==t.behind AND t.previous='" + chinese + "' AND e.pinyin='" + py + "' ORDER BY (e.probability + t.probability) DESC LIMIT 1;";
	//std::cout << sql << std::endl;
	const char* _sql = sql.c_str();
	if (sqlite3_prepare_v2(this->pDB, _sql, -1, &stmt, &zTail) == SQLITE_OK) {
		while (sqlite3_step(stmt) == SQLITE_ROW) {
			//std::vector<std::string> r;
			std::string _py(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0)));
			double pro = sqlite3_column_double(stmt, 1);
			//r.push_back(_py);
			//r.push_back(_prob);
			//res.push_back(r);
			res[_py] = pro;
		}
	}
	return res;
}

/*std::vector<std::vector<std::string>> HMMTable::QueryObservation() {

}*/
}