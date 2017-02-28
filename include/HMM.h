#pragma once
#include <cinttypes>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include "Query.h"
using namespace std;
//struct State2State
//{
//	int32_t state;
//	uint32_t cnt;
//	State2State(int32_t s, uint32_t c) :state(0), cnt(0) {};
//	bool operator() (const State2State &t) const
//	{
//		return state < t.state;
//	}
//};
struct Observed2State
{
	int32_t state;
	uint32_t cnt;
	Observed2State(int32_t s, uint32_t c) :state(0), cnt(0) {};
	bool operator() (const Observed2State &t) const
	{
		return state < t.state;
	}
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
	void insertMapRelations(string pinyin, wstring chars);
	void addCharacters(wstring characters,int cnt);
	void addChar2Char(wchar_t char1,wchar_t char2,int cnt);
	void addChar(wchar_t character,int cnt);
	size_t getPinyinIdFromCharId(size_t charId);
	size_t getPinyinIdFromPinyin(string pinyin);
	size_t getCharIdFromChar(wchar_t character);
	size_t getCharIdFromPinyinId(size_t pinyinId);
	void touchPinyinCharMap(size_t pinyinId, size_t charId);
	void addToPinyinCharMap(string pinyin,wstring chars);
	vector<Query> query(vector<string> pinyins);
private:
	vector<map<int32_t,uint32_t> > stateTransfers;
	vector<uint32_t> cntStateTransferOut;
	vector<vector<Observed2State> > observedTransfer;
	vector<int> charOccr; //single char occurence num
	vector<int> pinyinOccr; //single pinyin occurence num
	map<wchar_t, uint32_t> word2id;
	map<string, uint32_t> pinyin2id;
	vector<wchar_t> id2word;
	vector<string> id2pinyin;
	vector<set<uint32_t>> pinyinId2charsId;
	vector<uint32_t> charId2PinyinId;//ÉÐÎ´¿¼ÂÇ¶àÒô×Ö
private:
	vector<Query> Viterbi(vector<uint32_t> pinyinIds);
	uint32_t search(int deep, vector<uint32_t> &pyids,map<uint32_t,double> lastlayer,vector<uint32_t> &state);
};

