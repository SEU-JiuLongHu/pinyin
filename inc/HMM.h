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
#include "type.h"
using namespace std;
struct MatrixNode
{
	//stat: state id
	state stat;
	probability logp;
	//pre: record the previous state id which constitute the optimal path
	state pre;
	MatrixNode(state a = 0, probability b = numeric_limits<probability>::lowest()) :stat(a), logp(b) {};
};
typedef vector<MatrixNode> MatrixLayer;
typedef vector<MatrixLayer> Matrix;

struct HashFunc
{
	uint32_t operator()(const id_pair &key) const 
	{
		return hash<id>()(key.first)^hash<id>()(key.second);
	}
};
struct Collision
{
	bool operator() (const id_pair &lhs, const id_pair &rhs)
	{
		return (HashFunc()(lhs) == HashFunc()(rhs));
	}
};
struct Res
{
	probability prob;
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
	bool loadPinyinCharMap(string mapPath);
	void insert_PyChar_Relations(string pinyin, wstring chars);
	//void addCharacters(wstring characters,int cnt);
	//void addChar2Char(wchar_t char1,wchar_t char2,int cnt);
	//void addChar(wchar_t character,int cnt);
	//size_t getPinyinIdFromCharId(size_t charId);
	//size_t getPinyinIdFromPinyin(string pinyin);
	//size_t getCharIdFromChar(wchar_t character);
	//size_t getCharIdFromPinyinId(size_t pinyinId);
	//void touchPinyinCharMap(size_t pinyinId, size_t charId);
	//void addToPinyinCharMap(string pinyin,wstring chars);
	id get_charid(wchar_t ch);
	id get_pyid(string py);
	void add_char(wchar_t ch, uint32_t cnt);
	void add_chars(wstring chs, uint32_t cnt);
	vector<Res> query(vector<string> py,uint32_t topk);
	vector<Res> solve(uint32_t topk);
private:
	unordered_map<wchar_t, id> char2id;
	unordered_map<id, wchar_t> id2char;
	unordered_map<string, id> py2id;
	unordered_map<id, string> id2py;

	unordered_map<id,uint32_t> transf_out;
	unordered_map<id_pair, uint32_t,HashFunc> transf_to;
	unordered_map<id, uint32_t> char_freq;
	unordered_map<id, uint32_t> py_freq;
	unordered_map<id, id_list> pyid2charidlst;//map pingyinid into charid list 
	//ViterbiMatrix viterbi;
	Matrix matrix;
};

