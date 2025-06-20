// --- START OF FILE parmset.cpp (已修复) ---

#include "parmset.h"
#include "cexception.h"
#include <fstream>
#include <sstream>

using namespace std;

ParmSet::ParmSet() {
	//ctor
}

ParmSet::ParmSet(string parmfile) {
	string qq = parmfile;
	ifstream ips(parmfile.c_str());
	if (ips.fail())
		throw CException(std::string("Cannot open parameter file: ") + parmfile);
	ips >> *this;// input content in parmfile using ifstream ips. Ref. istream& operator>> (istream &is, ParmSet &ps)
}

ParmSet::~ParmSet() {
	//dtor
}

bool ParmSet::isDef(std::string key) {
	return table.find(key) != table.end();
}

bool ParmSet::getString(std::string kw, std::string &out) {
	map<string, string>::iterator itr = table.find(kw);
	if (itr == table.end()) return 0;
	out = itr->second;
	return 1;
}

bool ParmSet::getInt(std::string kw, int &out) {
	map<string, string>::iterator itr = table.find(kw);
	if (itr == table.end()) return 0;
	istringstream buf;
	buf.str(itr->second);
	buf >> out;
	return 1;
}

bool ParmSet::getFloat(std::string kw, double &out) {
	map<string, string>::iterator itr = table.find(kw);
	if (itr == table.end()) return 0;
	istringstream buf;
	buf.str(itr->second);
	buf >> out;
	return 1;
}

string ParmSet::getString(string key) {
	map<string, string>::iterator itr = table.find(key);
	if (itr == table.end()) return "";
	return itr->second;
}

int ParmSet::getInt(string key) {
	map<string, string>::iterator itr = table.find(key);
	if (itr == table.end()) return 0;
	istringstream buf;
	buf.str(itr->second);
	int v;
	buf >> v;
	// cout << "**ParmSet::getInt** " << key << "|" << v << endl; 
	return v;
}

double ParmSet::getFloat(string key) {
	map<string, string>::iterator itr = table.find(key);
	if (itr == table.end()) return 0;
	istringstream buf;
	buf.str(itr->second);
	double v;
	buf >> v;
	return v;
}

vector<double> ParmSet::getFloats(string key) {
	vector<double> vs;
	map<string, string>::iterator itr = table.find(key);
	if (itr == table.end()) return vs;
	string tmp = itr->second;
	size_t j;
	while ((j = tmp.find(',')) != string::npos) { tmp[j] = ' '; }
	//cerr << itr->second << endl;

	istringstream buf;
	buf.str(tmp);

	double v;
	while (buf >> v) vs.push_back(v);
	return vs;
}

ostream& operator<< (ostream &os, const ParmSet &ps) {
	map<string, string>::const_iterator itr = ps.table.begin();
	for (; itr != ps.table.end(); itr++) {
		os << itr->first << ":" << itr->second << endl;
	}
	return os;
}

// parse the input para file to std::map<std::string, std::string> table
// the key and value is seperated by '=', e.g., key = value, per line.
istream& operator>> (istream &is, ParmSet &ps) {
	const char* ws = " \f\n\r\t\v";
	string line;
	string s1;
	while (getline(is, line)) {
		string::size_type i, j, k;
		// # exist in line
		if ((j = line.find('#')) != string::npos) {
			if (j == 0) line.clear();// the first char is # --> comment line
			else line = line.substr(0, j - 1);
		}
		// cout << "**ParmSet**" << line << "    |   ";

		// remove blank char
		string::size_type pos = 0;
		while (pos < line.length()) {
			// 找到第一个非空白字符
			pos = line.find_first_not_of(" \t\n\v\f\r", pos);
			if (pos == string::npos) {
				break; // 如果没有找到非空白字符，退出循环
			}
			// 将非空白字符添加到结果字符串
			s1 += line[pos];
			// 移动到下一个字符
			++pos;
		}
		line = s1;// update line
		s1.clear();
		// cout << line << endl;

		// parse line to key & value by seperation of '='
		if ((i = line.find_first_not_of(ws)) != string::npos) {
			j = line.find_last_not_of(ws);
			if (j == string::npos) j = line.size();
			string token = line.substr(i, j - i + 1);
			k = token.find_first_of("=");
			if (k != string::npos) {
				string key, value;
				key = token.substr(0, k);
				value = token.substr(k + 1, token.size());
				// --- 修复后：去掉了这里多余的cout，它会干扰程序正常输出 ---
				// cout << key << "=" << value << endl;
				ps.table[key] = value;
			}
		}
	}
	return is;
} // --- 修复点：在这里添加了缺失的右大括号 ---

void ParmSet::putKeyValue(std::string kw, std::string vl) {
	table.insert(map<std::string, std::string>::value_type(kw, vl));
}

void ParmSet::changeKeyValue(std::string kw, std::string vl) {

	map<string, string>::iterator itr = table.find(kw);
	itr->second = vl;
}