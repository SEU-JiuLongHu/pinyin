#pragma once
#include "type.h"
#include <vector>
struct MatrixNode
{
	//stat: state id
	state stat;
	probability logp; 
	//pre: record the previous state id which constitute the optimal path
	state pre;
	MatrixNode(state a, probability b) :stat(a), logp(b) {};
};
typedef vector<MatrixNode> MatrixLayer;
typedef vector<MatrixLayer> Matrix;
class ViterbiMatrix
{
public:
	ViterbiMatrix();
	~ViterbiMatrix();
	void buildMatrix();
	void solve();
	void getTopK(size_t k);
private:
	Matrix matrix;
};

