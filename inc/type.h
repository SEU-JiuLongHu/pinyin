#pragma once
#include <cstdint>
#include <unordered_map>
#include <vector>
using namespace std;

namespace type {

	typedef uint32_t mhash;
	typedef uint32_t id;
	typedef uint32_t state;
	typedef uint32_t observation;
	typedef uint32_t count;
	typedef double probability;
	typedef vector<size_t> size_array;
	typedef pair<id, id> id_pair;
	typedef vector<id> id_list;


	struct HashIdPair
	{
		mhash operator()(const id_pair &key) const
		{
			return std::hash<id>()(key.first) ^ std::hash<id>()(key.second);
		}
	};
	struct Key
	{
		id pyid;
		id chid;
#ifdef DEBUG
		string pinyin;
		wchar_t ch;
		Key(id p,string py, id c,wchar_t cha) :pyid(p), chid(c),pinyin(py),ch(cha) {};
#endif // DEBUG
		Key(id p, id c) :pyid(p), chid(c) {};
		bool operator==(const Key &a) const
		{
			if (pyid == a.pyid&&chid == a.chid)
				return true;
			else
				return false;
		}
	};
	struct HashKey
	{
		mhash operator()(const Key &a) const
		{
			return hash<id>()(a.pyid) ^ hash<id>()(a.chid);
		}
	};
	typedef pair<Key, Key> Key_pair;
	struct HashKeyPair
	{
		mhash operator()(const Key_pair &a) const
		{
			return HashKey()(a.first) ^ HashKey()(a.second);
		}
	};



	struct MatrixNode
	{
#ifdef DEBUG
		wchar_t ch;
#endif // DEBUG
		//stat: state id
		state stat;
		probability logp;
		//pre: record the previous state id which constitute the optimal path
		state pre;
#ifdef DEBUG
		MatrixNode(wchar_t c, state a = 0, probability b = numeric_limits<probability>::lowest()) :ch(c),stat(a), logp(b) {};
#endif // DEBUG
		MatrixNode(state a = 0, probability b = numeric_limits<probability>::lowest()) :stat(a), logp(b) {};

	};
	typedef vector<MatrixNode> MatrixLayer;
	typedef vector<MatrixLayer> Matrix;
}