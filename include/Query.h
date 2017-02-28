#pragma once
#include <string>
using namespace std;
struct Query
{
	wstring result;
	double probability;
	bool operator() (const Query &query) const
	{
		return probability > query.probability;
	}
};