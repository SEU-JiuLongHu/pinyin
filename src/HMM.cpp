#include "HMM.h"
#include "PinyinConverter.h"
//#define _check_vector_size(x) if (x.size()>x.capacity()*0.5) x.reserve(2*x.size());
HMM::HMM()
{
#ifdef DEBUG

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
		//std::ios::sync_with_stdio(false);
		string wordUtf8,pinyin;
		wstring wordUnicode;
		uint32_t cnt;
		while (fin >> wordUtf8)
		{
			wordUnicode = PinyinConverter::Utf8ToUnicode(wordUtf8);
			for (size_t i = 0; i < wordUnicode.size(); i++)
			{
				fin >> pinyin;
			}
			fin >> cnt;
			if (wordUnicode.size() > 1)
				add_chars(wordUnicode, cnt);
			else if (wordUnicode.size() == 1)
				add_char(wordUnicode[0], cnt);
			else continue;
		}
		return true;
	}
	else
		return false;
}
id HMM::get_charid(wchar_t ch)
{
	auto iter= char2id.find(ch);
	if (iter == char2id.end())
	{
		id charid = char2id.size();
		char2id[ch] = charid;
		id2char[charid] = ch;
		return charid;
	}
	else
		return (*iter).second;
}
id HMM::get_pyid(string py)
{
	auto iter = py2id.find(py);
	if (iter == py2id.end())
	{
		id pyid = py2id.size();
		py2id[py] = pyid;
		id2py[pyid] = py;
		return pyid;
	}
	else
		return (*iter).second;
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
			insert_PyChar_Relations(pinyin, chars);
		}
		return true;
	}
	else
		return false;
}
//建立拼音id与字id之间的映射关系
void HMM::insert_PyChar_Relations(string pinyin, wstring chars)
{
	id pyid=get_pyid(pinyin);
	if (pyid2charidlst.find(pyid) == pyid2charidlst.end())
		pyid2charidlst[pyid] = id_list();
	for (auto &ch : chars )
	{
		id charid = get_charid(ch);
		pyid2charidlst[pyid].push_back(charid);
	}

}

void HMM::add_char(wchar_t ch, uint32_t cnt)
{
	id i = get_charid(ch);
	char_freq[i] += cnt;
}

void HMM::add_chars(wstring chs, uint32_t cnt)
{
	wchar_t pre;
	for (auto iter=chs.begin();iter!=chs.end();iter++)
	{
		if (iter != chs.begin())
		{
			id i = get_charid(pre);
			id ii = get_charid(*iter);
			id_pair pair = make_pair(i, ii);
			transf_to[pair] += cnt;
			transf_out[i] += cnt;
		}
		pre = *iter;
		add_char(pre,cnt);
	}
}

vector<Res> HMM::query(vector<string> pinyins,uint32_t topk)
{
	matrix.clear();
	matrix.resize(pinyins.size());
	for (size_t i=0;i<pinyins.size();i++)
	{
		string py = pinyins[i];
		id pyid = py2id[py];
		auto charidlst = pyid2charidlst[pyid];
		size_t charsize = charidlst.size();
		matrix[i].resize(charsize);
		/*
		没有多音字，暂时使用字符的频率和计算第一个字的概率
		*/
		int sum = 0;
		for (size_t j = 0; j < charsize; j++)
		{
			state curstat = charidlst[j];
			sum += char_freq[curstat];
		}
		for (size_t j = 0; j < charsize; j++)
		{
			state curstat = charidlst[j];
			matrix[i][j] = MatrixNode(curstat);
			if (!i)
			{
				matrix[i][j].logp = log(1.0*(char_freq[curstat]+1) / (sum+charsize));
			}
			else
			{
				state prestat;
				probability prob;
				for (size_t k = 0; k < matrix[i - 1].size(); k++)
				{
					prestat = matrix[i - 1][k].stat;
					prob = matrix[i - 1][k].logp + log(1.0*(transf_to[make_pair(prestat, curstat)]+1) / (transf_out[prestat]+ char2id.size()));
					if (prob > matrix[i][j].logp)
					{
						matrix[i][j].logp = prob;
						matrix[i][j].pre = k;
					}
				}
			}
		}
	}
	return solve(topk);
}
vector<Res> HMM::solve(uint32_t topk)
{
	static auto cmp = [](const MatrixNode &a, const MatrixNode &b) {
		return a.logp > b.logp;
	};
	assert(matrix.size()>0);
	uint32_t L = matrix.size() - 1;
	sort(matrix[L].begin(), matrix[L].end(), cmp);
	vector<Res> records;
	for (int i = 0; i < min(topk,matrix[L].size()); i++)
	{
		int cur=i;
		Res record(matrix[L][i].logp);
		for (int j = L; j >=0 ; j--)
		{
			record.str.push_back(id2char[matrix[j][cur].stat]);
			cur = matrix[j][cur].pre;
		}
		reverse(record.str.begin(), record.str.end());
		records.push_back(record);
	}
	return records;
}