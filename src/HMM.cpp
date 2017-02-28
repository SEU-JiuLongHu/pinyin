#include "HMM.h"
#include "PinyinConverter.h"
#define _check_vector_size(x) if (x.size()>x.capacity()*0.5) x.reserve(2*x.size());
HMM::HMM()
{
#ifdef DEBUG
	stateTransfers.reserve(30000);
	observedTransfer.reserve(1000);
	charOccr.reserve(30000);
	pinyinOccr.reserve(1000);
	id2word.reserve(30000);
	id2pinyin.reserve(1000);
	pinyinId2charsId.reserve(1000);
	charId2PinyinId.reserve(1000);
	cntStateTransferOut.reserve(30000);
#endif // DEBUG
}

HMM::HMM(string corpus)
{
	


}
HMM::~HMM()
{
}

bool HMM::loadCorpus(string corpus)
{
	ifstream fin(corpus);
	if (fin)
	{
		string wordUtf8;
		wstring wordUnicode;
		int cnt;
		while (fin >> wordUtf8 >> cnt)
		{

			wordUnicode = PinyinConverter::Utf8ToUnicode(wordUtf8);
			
			if (wordUnicode.size() > 1)
				addCharacters(wordUnicode, cnt);
			else if (wordUnicode.size() == 1)
				addChar(wordUnicode[0], cnt);
			else continue;
		}
		return true;
	}
	else
		return false;
}
bool HMM::loadPinyinCharMap(string mapPath)
{
	ifstream fin(mapPath);
	if (fin)
	{
		string pinyin;
		string wordsUtf8;
		wstring chars;
		while (fin >> pinyin >> wordsUtf8)
		{
			chars = PinyinConverter::Utf8ToUnicode(wordsUtf8);
			insertMapRelations(pinyin, chars);
		}
		return true;
	}
	else
		return false;
}
//建立拼音id与字id之间的映射关系
void HMM::insertMapRelations(string pinyin, wstring chars)
{
	size_t pinyinId, charId;
	pinyinId = getPinyinIdFromPinyin(pinyin);
	_check_vector_size(pinyinId2charsId);
	while (pinyinId >= pinyinId2charsId.size())
		pinyinId2charsId.push_back(set<uint32_t>());
	for (size_t i = 0; i < chars.size(); i++)
	{
		charId = getCharIdFromChar(chars[i]);
		_check_vector_size(charId2PinyinId);
		while (charId >= charId2PinyinId.size()) charId2PinyinId.push_back(0);
		charId2PinyinId[charId] = pinyinId;
		_check_vector_size(pinyinId2charsId);
		while (pinyinId >= pinyinId2charsId.size()) pinyinId2charsId.push_back(set<uint32_t>());
		pinyinId2charsId[pinyinId].insert(charId);
	}
}
void HMM::addChar2Char(wchar_t char1, wchar_t char2, int cnt)
{
	size_t id1 = getCharIdFromChar(char1),id2= getCharIdFromChar(char2);
	stateTransfers[id1][id2] += cnt;
	cntStateTransferOut[id1] += cnt;
	//
	//addChar(char1, cnt);
	//addChar(char2, cnt);
}
void HMM::addChar(wchar_t character, int cnt)
{
	size_t id = getCharIdFromChar(character);
	charOccr[id] += cnt;
	size_t pyId = getPinyinIdFromCharId(id);
	pinyinOccr[pyId] += cnt;
}

size_t HMM::getPinyinIdFromCharId(size_t charId)
{
	assert(charId<charId2PinyinId.size());
	return charId2PinyinId[charId];
}

size_t HMM::getPinyinIdFromPinyin(string pinyin)
{
	if (pinyin2id.find(pinyin) == pinyin2id.end())
	{
		pinyin2id[pinyin] = pinyin2id.size();
		_check_vector_size(id2pinyin);
		id2pinyin.push_back(pinyin);
		_check_vector_size(pinyinOccr);
		while (pinyinOccr.size() < pinyin2id.size())
			pinyinOccr.push_back(0);
		return pinyin2id.size() - 1;
	}
	else
		return pinyin2id[pinyin];
}

size_t HMM::getCharIdFromChar(wchar_t character)
{
	if (word2id.find(character) == word2id.end())
	{
		word2id[character]= word2id.size();
		id2word.push_back(character);
		if (stateTransfers.size() < word2id.size())
		{
			stateTransfers.push_back(map<int32_t, uint32_t>());
			cntStateTransferOut.push_back(0);
		}
		if (charOccr.size() < word2id.size())
			charOccr.push_back(0);
		return word2id.size() - 1;
	}
	else
		return word2id[character];
}

size_t HMM::getCharIdFromPinyinId(size_t pinyinId)
{
	return size_t();
}

void HMM::touchPinyinCharMap(size_t pinyinId, size_t charId)
{
	//pinyinId2charsId.size();


	//
}

void HMM::addToPinyinCharMap(string pinyin, wstring chars)
{
	//int pinyinId = getPinyinIdFromPinyin(pinyin);
	//int charId;
	//for (size_t i = 0; i < chars.size(); i++)
	//{
	//	charId = getCharIdFromChar(chars[i]);
	//	touchPinyinCharMap(pinyinId, charId);
	//}
}

vector<Query> HMM::query(vector<string> pinyins)
{
	vector<uint32_t> ids;
	for (size_t i = 0; i < pinyins.size(); i++)
	{
		if (pinyin2id.find(pinyins[i]) == pinyin2id.end())
			throw string("不含该拼音");
		else
			ids.push_back(pinyin2id[pinyins[i]]);
	}

	return Viterbi(ids);
}
vector<Query> HMM::Viterbi(vector<uint32_t> pinyinIds)
{
	double maxProb=0,currProb=1;
	vector<uint32_t> states(pinyinIds.size());
	search(0, pinyinIds,map<uint32_t,double>(),states);
	wstring ans;
	for (size_t i = 0; i < states.size(); i++)
		ans+=id2word[states[i]];
	ofstream fo("out.txt");
	fo << PinyinConverter::UnicodeToUtf8(ans);
	fo.close();
	return vector<Query>();
}
uint32_t HMM::search(int deep,vector<uint32_t> &pyids,map<uint32_t,double> lastlayer,vector<uint32_t> &result)
{
	if (deep >= pyids.size())
	{
		uint32_t maxProbState;
		double p = 0;
		for (map<uint32_t,double>::iterator i=lastlayer.begin();i!=lastlayer.end();i++)
		{
			if ((*i).second > p)
			{
				maxProbState = (*i).first;
				p = (*i).second;
			}
		}
		return maxProbState;
	}
	//auto max = [](const double &a, const double &b) {return a > b; };
	map<uint32_t, double> currlayer;
	if (!deep)
	{
		double pstate2observed;
		for (set<uint32_t>::iterator i = pinyinId2charsId[pyids[deep]].begin(); i != pinyinId2charsId[pyids[deep]].end(); i++)
		{
				pstate2observed =double (charOccr[*i] + 1) / (pinyinOccr[pyids[deep]] + pinyinId2charsId[pyids[deep]].size());
				currlayer[*i] = pstate2observed;
		}
		int maxProbstate=search(deep + 1, pyids,currlayer, result);
		result[deep] = maxProbstate;
		return 0;
	}
	else
	{
		double pstate2observed, pstate2state;
		uint32_t currstate, prestate;
		map<uint32_t, uint32_t> maxProbPre;
		for (set<uint32_t>::iterator i = pinyinId2charsId[pyids[deep]].begin(); i != pinyinId2charsId[pyids[deep]].end(); i++)
		{
			currstate = *i;
			currlayer[currstate] = 0;
			for (map<uint32_t, double>::iterator j = lastlayer.begin(); j != lastlayer.end(); j++)
			{
				prestate = (*j).first;
				//pstate2observed =double (charOccr[*i] + 1) / (pinyinOccr[pyids[deep]] + pinyinId2charsId[pyids[deep]].size());
				pstate2observed = 1;
				pstate2state = double(stateTransfers[prestate][currstate] + 1) / (cntStateTransferOut[prestate] + stateTransfers[prestate].size());
				if (lastlayer[prestate] * pstate2state*pstate2observed > currlayer[currstate])
				{
					currlayer[currstate] = lastlayer[prestate] * pstate2state*pstate2observed;
					maxProbPre[currstate] = prestate;
				}
			}
		}
		int maxProbstate=search(deep + 1, pyids,currlayer, result);
		result[deep] = maxProbstate;
		return maxProbPre[maxProbstate];
	}

}
void HMM::addCharacters(wstring word, int cnt)
{
	for (int i = 0; i < word.size() - 1; i++)
	{
		addChar2Char(word[i], word[i + 1], cnt);
	}
}
