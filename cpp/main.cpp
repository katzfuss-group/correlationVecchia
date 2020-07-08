#include "SortSparse.h"

int main()
{
	Node A; A.id = 0; A.rank = 2; A.val = 1.5;
	Node B; B.id = 1; B.rank = 1; B.val = 0.5;
	Node C; C.id = 2; C.rank = 0; C.val = 2.5;
	Node D; D.id = 3; D.rank = 3; D.val = 5.5;
	Node E; E.id = 4; E.rank = 4; E.val = 6.5;
	Node F; F.id = 5; F.rank = 5; F.val = 7.5;

	MutHeap H; H.nodes = { A, B, C }; H.lookup = { 0, 1, 2 };

	for (int i = 0; i < 3; i++) {
		cout << "value = " << double(H.nodes[i].val) << '\n';
		cout << "lookup = " << int(H.lookup[i]) << '\n';
		cout << '\n';
	}

	_swap(&H, 1, 2);

	for (int i = 0; i < 3; i++) {
		cout << "value = " << double(H.nodes[i].val) << '\n';
		cout << "lookup = " << int(H.lookup[i]) << '\n';
		cout << '\n';
	}

	cout << '\n' << '\n' << '\n' << '\n';

	for (int i = 0; i < 3; i++) {
		cout << "value = " << double(H.nodes[i].val) << '\n';
		cout << "lookup = " << int(H.lookup[i]) << '\n';
		cout << '\n';
	}

	signed int a = _moveDown(&H, 0);
	cout << signed int(a) << '\n' << '\n';

	for (int i = 0; i < 3; i++) {
		cout << "value = " << double(H.nodes[i].val) << '\n';
		cout << "lookup = " << int(H.lookup[i]) << '\n';
		cout << '\n';
	}

	cout << '\n' << '\n' << '\n' << '\n';

	Node tempNode = topNode(&H);
	cout << "id = " << signed int(tempNode.id) << '\n';
	cout << "val = " << double(tempNode.val) << '\n';
	cout << "rank = " << signed int(tempNode.rank) << '\n';

	tempNode = topNode_rankUpdate(&H);
	cout << "id = " << signed int(tempNode.id) << '\n';
	cout << "val = " << double(tempNode.val) << '\n';
	cout << "rank = " << signed int(tempNode.rank) << '\n';

	cout << '\n' << '\n' << '\n' << '\n';

	double b = update(&H, 0, 0);
	tempNode = topNode(&H);
	for (int i = 0; i < 3; i++) {
		cout << "value = " << double(H.nodes[i].val) << '\n';
		cout << "lookup = " << int(H.lookup[i]) << '\n';
		cout << '\n';
	}

	cout << '\n' << '\n' << '\n' << '\n';

	H.nodes.insert(H.nodes.end(), D);
	H.lookup.insert(H.lookup.end(), 3);
	for (signed int i = 0; i < H.nodes.size(); i++) {
		cout << "value = " << double(H.nodes[i].val) << '\n';
		cout << "lookup = " << int(H.lookup[i]) << '\n';
		cout << '\n';
	}

	cout << '\n' << '\n' << '\n' << '\n';

	vector<Node> G = { E, F };
	vector<signed int> lookdown = { 4, 5 };

	// Node G[2] = { E, F };
	// signed int lookdown[2] = { 4, 5 };

	H.nodes.insert(H.nodes.end(), G.begin(), G.end());
	H.lookup.insert(H.lookup.end(), lookdown.begin(), lookdown.end());
	for (signed int i = 0; i < H.nodes.size(); i++) {
		cout << "value = " << double(H.nodes[i].val) << '\n';
		cout << "lookup = " << int(H.lookup[i]) << '\n';
		cout << '\n';
	}

	cout << '\n' << '\n' << '\n' << '\n' << "Member";

	Member mA = assignMember(9.5, 0);
	Member mB = assignMember(11.5, 1);

	cout << '\n' << bool(mA < mB);

	cout << '\n' << '\n' << '\n' << '\n';

	vector<int> tempVec;

	tempVec.resize(3, 100);
	// tempVec[0] = 1;
	// tempVec[1] = 2;
	// tempVec[2] = 3;

	for (signed int i = 0; i < tempVec.size(); i++) {
		cout << " " << int(tempVec[i]);
	}

	cout << '\n' << '\n' << '\n' << '\n';

	vector<Member> tempVecMem;
	tempVecMem.resize(10);
	for (signed int i = 0; i < 10; i++) {
		tempVecMem.at(i).id = i;
		tempVecMem.at(i).val = 0.1 * i;
	}

	tempVecMem = subMember(tempVecMem, 0, 7);
	for (signed int i = 0; i < tempVecMem.size(); i++) {
		cout << signed int(tempVecMem[i].id) << '\n';
	}

	cout << '\n' << '\n' << '\n' << '\n';

	// Data X;
	// vector<double> row1 = { 0.768447675, 0.673958695, 0.313243956, 0.586022124, 0.268639569, 0.163665819, 0.865412143, 0.2856979, 0.275819115, 0.58231778, 0.705860332, 0.281065895, 0.209230168, 0.614254771, 0.555668301, 0.479999833, 0.356221085, 0.529253324, 0.900681479, 0.621378715, 0.570613189, 0.37498027, 0.1917802, 0.097669803, 0.946696535, 0.122594436, 0.514679883, 0.03347555, 0.301603402, 0.458475706 };
	// vector<double> row2 = { 0.940515001, 0.395453112, 0.662554816, 0.052133163, 0.108870741, 0.473016816, 0.617491888, 0.463847208, 0.446568065, 0.255981303, 0.291978267, 0.792931029, 0.918165131, 0.802664694, 0.94078227, 0.790200953, 0.900924596, 0.031830968, 0.940299242, 0.348172765, 0.203996623, 0.759754654, 0.234543676, 0.627092997, 0.54690429, 0.110364564, 0.628172302, 0.707863993, 0.31084741, 0.46721859 };

	// X.entries.push_back(row1);
	// X.entries.push_back(row2);
	// X.colSize = 30; X.rowSize = 2;

	cout << double(x.entries[0][0]);

	cout << '\n' << '\n' << '\n' << '\n';

	cout << double(dist2Func(1, 1)) << '\n' << double(dist2Func(1, 2));

	cout << '\n' << '\n' << '\n' << '\n';

	output result = sortSparse(x.colSize, 3.5, dist2Func, 0);

	for (signed int i = 0; i < result.P.size(); i++) {
		cout << '\n' << signed int(result.P[i]);
	}

	cout << '\n' << '\n';

	for (signed int i = 0; i < result.revP.size(); i++) {
		cout << '\n' << signed int(result.revP[i]);
	}

	cout << '\n' << '\n';

	for (signed int i = 0; i < result.colptr.size(); i++) {
		cout << '\n' << signed int(result.colptr[i]);
	}

	cout << '\n' << '\n';

	for (signed int i = 0; i < result.distances.size(); i++) {
		cout << '\n' << double(result.distances[i]);
	}

	/*
	cout << '\n' << '\n' << '\n' << '\n';

	for (signed int i = 0; i < result.checkInt.size(); i++) {
		cout << '\n' << signed int(result.checkInt[i]);
	}

	cout << '\n' << '\n';

	for (signed int i = 0; i < result.checkDouble.size(); i++) {
		cout << '\n' << double(result.checkDouble[i]);
	}
	*/

	return(0);
}

// 0.768447675, 0.673958695, 0.313243956, 0.586022124, 0.268639569, 0.163665819, 0.865412143, 0.2856979, 0.275819115, 0.58231778, 0.705860332, 0.281065895, 0.209230168, 0.614254771, 0.555668301, 0.479999833, 0.356221085, 0.529253324, 0.900681479, 0.621378715, 0.570613189, 0.37498027, 0.1917802, 0.097669803, 0.946696535, 0.122594436, 0.514679883, 0.03347555, 0.301603402, 0.458475706 

// 0.940515001, 0.395453112, 0.662554816, 0.052133163, 0.108870741, 0.473016816, 0.617491888, 0.463847208, 0.446568065, 0.255981303, 0.291978267, 0.792931029, 0.918165131, 0.802664694, 0.94078227, 0.790200953, 0.900924596, 0.031830968, 0.940299242, 0.348172765, 0.203996623, 0.759754654, 0.234543676, 0.627092997, 0.54690429, 0.110364564, 0.628172302, 0.707863993, 0.31084741, 0.46721859