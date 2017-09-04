#pragma once
#include <cinttypes>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <sstream>
#include "type.h"
using namespace std;
using namespace type;

struct Res
{
	type::probability prob;
	wstring str;
	Res(probability p) :prob(p) {};
};
class HMM
{
public:
	HMM();
	HMM(string corpus);
	~HMM();
	//vector tokenizer(vector<wstring> input);
	bool loadCorpus(string corpusPath);
	
	id get_charid(wchar_t ch);
	id get_pyid(string py);
	void add_char(wchar_t ch, string pinyin, type::count cnt);
	void add_chars(wchar_t character1, string pinyin1, wchar_t character2, string pinyin2, type::count cnt);
	vector<Res> query(vector<string> py,uint32_t topk);
	vector<Res> solve(uint32_t topk);
private:
	void init_pinyin2chars_table();
	void insert_pychar_relations(string pinyin, wstring chars);
private:
	unordered_map<wchar_t, id> char2id;
	unordered_map<id, wchar_t> id2char;
	unordered_map<string, id> py2id;
	unordered_map<id, string> id2py;

	//用来计算转移概率
	unordered_map<Key, type::count, HashKey> transf_out;
	unordered_map<Key_pair, type::count, HashKeyPair> transf_to; 

	//用来计算初始概率
	unordered_map <Key, type::count, HashKey > char_freq;//the count of a char given the pinyin
	unordered_map <id, type::count> py_freq;//the count of the pinyin

	//用来计算发射概率
	unordered_map<id_pair, type::count, HashIdPair> emit_to;//从字转移到拼音的数量
	unordered_map<id, type::count> emit_out;//字的数量
	
	unordered_map<id, id_list> pyid2charidlst;//map pingyinid into charid list 
	//ViterbiMatrix viterbi;
	Matrix matrix;

	static const map<string, wstring> pinyin2chars;
};
