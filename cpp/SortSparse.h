#pragma once

#include "Matrix.h"
#include "MutHeap.h"
#include "math.h"
#include "algorithm"

struct Member
{
	double val;
	signed int id;

	bool operator<(const Member& other);
};

bool Member::operator<(const Member& other)
{
	/*
	if (id == other.id && val == other.val) {
		return(false);
	}
	else if (id >= other.id && val <= other.val) {
		return(true);
	}
	else {
		return(false);
	}
	*/

	if (val < other.val && id > other.id ) {
		return(true);
	}
	else {
		return(false);
	}
}

bool compareMember(Member const a, Member const b)
{
	/*
	if (a.id == b.id && a.val == b.val) {
		return(false);
	}
	else if (a.id >= b.id && a.val <= b.val) {
		return(true);
	}
	else {
		return(false);
	}
	*/

	if (a.val < b.val) { // && a.id > b.id
		return(true);
	}
	else {
		return(false);
	}
}

Member assignMember(double val, signed int id)
{
	Member output = { val, id };
	return(output);
}

/*
bool islessMember(Member const a, Member const b)
{
	if (a.id == b.id && a.val == b.val) {
		return(false);
	}
	else if (a.id >= b.id && a.val <= b.val) {
		return(true);
	}
	else {
		return(false);
	}
}
*/

/*
bool operator<(Member const a, Member const b)
{
	if (a.id < b.id && a.val < b.val) {
		return(true);
	}
	else {
		return(false);
	}
}
*/

/*
bool operator <(Member const a, Member const b)
{
	if (a.id == b.id && a.val == b.val) {
		return(false);
	}
	else if (a.id <= b.id && a.val <= b.val) {
		return(true);
	}
	else {
		return(false);
	}
}
*/

struct ChildList
{
	signed int NParents = 0, NChildren = 0, NBuffer = 0;

	vector<signed int> P; // This array gives for contains the ordering.The i-th parent in the daycare has id P[i]
	vector<signed int> revP; // This array contains as the i-th element the number that the ith parent has with respect to the multiresolution ordering.
	vector<signed int> colptr; // The array that contains the first "child" for every parent
	vector<Member> rowval; // The array that contains the global id-s of the children
};

void newParent(ChildList *dc, signed int idParent)
{
	dc->P.at(dc->NParents) = idParent;
	dc->revP.at(idParent) = dc->NParents;

	dc->NParents++;

	dc->colptr.at(dc->NParents - 1) = dc->NChildren;
	dc->colptr.at(dc->NParents) = dc->NChildren;
}

void newChild(ChildList *dc, Member child)
{
	if (dc->NChildren >= dc->NBuffer) {
		if (dc->NChildren <= (signed int)(1e6)) {
			dc->NBuffer = 2 * dc->NBuffer;
		}
		else {
			dc->NBuffer = dc->NBuffer + (signed int)(1e6);
		}
	}

	dc->rowval.reserve(dc->NBuffer);

	dc->colptr.at(dc->NParents)++;
	dc->rowval.at(dc->NChildren) = child;

	dc->NChildren++;
	// dc->rowval.insert(dc->rowval.end(), child);
}

void newChildren(ChildList *dc, vector<Member> children)
{
	while (dc->NChildren + children.size() >= dc->NBuffer - 1) {

		if (dc->NChildren <= (signed int)(1e6)) {
			dc->NBuffer = 2 * dc->NBuffer;
		}
		else {
			dc->NBuffer = dc->NBuffer + (signed int)(1e6);
		}
	}

	dc->rowval.reserve(dc->NBuffer);

	dc->NChildren += children.size();
	dc->colptr.at(dc->NParents) += children.size();
		
	dc->rowval.resize(dc->NChildren);
	for (signed int i = dc->NChildren - children.size(); i < dc->NChildren; i++) {
		dc->rowval.at(i) = children.at(i - dc->NChildren + children.size());
	}	

	// dc->rowval.insert(dc->rowval.end(), children.begin(), children.end());
}

vector<Member> subMember(vector<Member> vec, signed int a, signed int b)
{
	auto first = vec.begin() + a;
	auto last = vec.begin() + b;

	vector<Member> subvec(first, last);

	return(subvec);
}

void _determineChildren(MutHeap *h, ChildList *dc, vector<Member> *parents, Node pivot, vector<Member> *buffer, double rho, double dist2Func(signed int, signed int))
{
	double distToParent = parents->at(pivot.id).val;
	double lengthScale = pivot.val;
	signed int iterBuffer = 0;

	Member candidate;
	double dist, dist2, newDist;
	vector<Member> viewBuffer;	

	signed int start = dc->colptr[dc->revP[parents->at(pivot.id).id]];
	signed int end = dc->colptr[dc->revP[parents->at(pivot.id).id] + 1];

	// printf("first i = %10d\n", start); printf("last i = %10d\n", end);

	for (signed int i = start; i < end; i++) {

		candidate = dc->rowval.at(i);	
		dist2 = dist2Func(candidate.id, pivot.id); // printf("\n\ndist2 = %10f\n", dist2); printf("%10d\n", i); printf("candidate.id = %10d /", candidate.id); printf("pivot.id = %10d \n", pivot.id);

		if (dc->revP.at(candidate.id) == -1 && dist2 <= pow(lengthScale * rho, 2.0)) {

			dist = sqrt(dist2); // printf("dist = %10f\n", dist);

			buffer->at(iterBuffer) = assignMember(dist, candidate.id);
			iterBuffer++;

			newDist = update(h, candidate.id, dist);

			// printf("newDist = %10f\n", newDist); printf("h0 = %10f\n", h->nodes[0].val);

			if (dist + rho * newDist <= rho * lengthScale && dist < parents->at(candidate.id).val) {
				parents->at(candidate.id) = assignMember(dist, pivot.id);
			}
			
		}
		else if (candidate.val > distToParent + lengthScale * rho) {
			break;
		}
	}

	/*
	for (signed int i = dc->colptr[dc->revP[parents->at(pivot.id).id]]; i < dc->colptr[dc->revP[parents->at(pivot.id).id] + 1]; i++) {

		printf("first i = %10d\n", dc->colptr[dc->revP[parents->at(pivot.id).id]]); printf("last i = %10d\n", dc->colptr[dc->revP[parents->at(pivot.id).id] + 1]);

		candidate = dc->rowval.at(i);	
		dist2 = dist2Func(candidate.id, pivot.id); printf("\n\ndist2 = %10f\n", dist2); printf("%10d\n", i); printf("candidate.id = %10d /", candidate.id); printf("pivot.id = %10d \n", pivot.id);

		if (dc->revP.at(candidate.id) == -1 && dist2 <= pow(lengthScale * rho, 2.0)) {

			dist = sqrt(dist2); printf("dist = %10f\n", dist);

			buffer->at(iterBuffer) = assignMember(dist, candidate.id);
			iterBuffer++;

			newDist = update(h, candidate.id, dist);

			printf("newDist = %10f\n", newDist); printf("h0 = %10f\n", h->nodes[0].val);
			// break;

			if (dist + rho * newDist <= rho * lengthScale && dist < parents->at(candidate.id).val) {
				parents->at(candidate.id) = assignMember(dist, pivot.id);
			}
			else if (candidate.val > distToParent + lengthScale * rho) {
				break;
			}
		}
	}
	*/

	/*
	for (signed int i = dc->colptr[dc->revP[parents->at(pivot.id).id]]; i < dc->colptr[dc->revP[parents->at(pivot.id).id] + 1]; i++) {

		candidate = dc->rowval.at(i);	
		dist2 = dist2Func(candidate.id, pivot.id); printf("\n\ndist2 = %10f\n", dist2);

		if (dc->revP.at(candidate.id) == -1 && dist2 <= pow(lengthScale * rho, 2.0)) {

			dist = sqrt(dist2); printf("dist = %10f\n", dist);

			buffer->at(iterBuffer) = assignMember(dist, candidate.id);
			iterBuffer++;

			newDist = update(h, candidate.id, dist);

			printf("newDist = %10f\n", newDist); printf("h0 = %10f\n", h->nodes[0].val);
			// break;

			if (dist + rho * newDist <= rho * lengthScale && dist < parents->at(candidate.id).val) {
				parents->at(candidate.id) = assignMember(dist, pivot.id);
			}

		}
		else if (candidate.val > distToParent + lengthScale * rho) {

			break;
		}

	}
	*/

	/*
	for (signed int i = dc->colptr[dc->revP[parents->at(pivot.id).id]]; i < dc->colptr[dc->revP[parents->at(pivot.id).id] + 1]; i++) {

		candidate = dc->rowval.at(i);
		dist2 = dist2Func(candidate.id, pivot.id);

		if (dc->revP.at(candidate.id) == -1 && dist2 <= pow(lengthScale * rho, 2.0)) {

			dist = sqrt(dist2);

			buffer->at(iterBuffer) = assignMember(dist, candidate.id);
			newDist = update(h, candidate.id, dist);

			iterBuffer++;

			if (dist + rho * newDist <= rho * lengthScale && dist < parents->at(candidate.id).val) {
				parents->at(candidate.id) = assignMember(dist, pivot.id);
				
			}
			else if (candidate.val > distToParent + lengthScale * rho) {
				break;
			}
		}
	}
	*/

	viewBuffer = subMember(*buffer, 0, iterBuffer); // printf("%10d", iterBuffer); printf("%10d", viewBuffer.size());
	sort(viewBuffer.begin(), viewBuffer.end(), compareMember);
	newParent(dc, pivot.id); printf("%10d", pivot.id);
	newChildren(dc, viewBuffer);
}	

struct output
{
	vector<signed int> colptr;
	vector<Member> rowval;
	vector<signed int> P;
	vector<signed int> revP;
	vector<double> distances;

	// vector<signed int> checkInt;
	// vector<double> checkDouble;
};

output sortSparse(signed int N, double rho, double dist2Func(signed int, signed int), signed int initInd)
{
	MutHeap h;
	ChildList dc;
	vector<Member> nodeBuffer;
	vector<double> distances;

	vector<Member> viewBuffer;
	vector<Member> parents;

	output result;

	h.nodes.resize(N);
	h.lookup.resize(N);

	for (signed int i = 0; i < N; i++) {
		h.nodes[i].val = numeric_limits<double>::max();
		h.nodes[i].id = i;
		h.nodes[i].rank = 0;
		h.lookup[i] = i;
	}

	dc.NParents = 0; dc.NChildren = 0; dc.NBuffer = N;
	dc.P.resize(N, -1); dc.revP.resize(N, -1); dc.colptr.resize(N + 1, -1);
	dc.rowval.resize(N);

	nodeBuffer.resize(N);
	distances.resize(N, -1.0);
	
	newParent(&dc, initInd);
	h.nodes.at(initInd).rank = numeric_limits<signed int>::max();
	distances.at(0) = numeric_limits<double>::max();
	
	for (signed int i = 0; i < N; i++) {
		nodeBuffer[i].val = update(&h, i, sqrt(dist2Func(i, initInd)));
		nodeBuffer[i].id = i;
	}

	/*
	for (signed int j = 0; j < h.nodes.size(); j++) {
		printf("hnodes.id  = %10d\n", h.nodes.at(j).id);
		printf("hnodes.val = %10f\n", h.nodes.at(j).val);
	}

	for (signed int j = 0; j < N; j++) {
		_moveDown(&h, j);
	}

	printf("\n========================================================================\n");

	for (signed int j = 0; j < h.nodes.size(); j++) {
		printf("hnodes.id  = %10d\n", h.nodes.at(j).id);
		printf("hnodes.val = %10f\n", h.nodes.at(j).val);
	}
	*/

	viewBuffer = subMember(nodeBuffer, 0, N);
	sort(viewBuffer.begin(), viewBuffer.end(), compareMember);
	newChildren(&dc, viewBuffer);
	
	parents.resize(N);
	for (signed int i = 0; i < N; i++) {
		parents[i] = assignMember(sqrt(dist2Func(initInd, i)), initInd);
	}

	for (signed int i = 1; i < N; i++) {
		distances[i] = topNode_rankUpdate(&h).val; 
		_determineChildren(&h, &dc, &parents, topNode(&h), &nodeBuffer, rho, dist2Func);

		/*
		if (i == 1) {
			printf("topNode.id = %10d\n", topNode(&h).id);
			for (signed int j = 0; j < dc.rowval.size(); j++) {
				printf("rowval.id  = %10d\n", dc.rowval.at(j).id);
				printf("rowval.val = %10f\n", dc.rowval.at(j).val);
			}
			printf("\n");
			for (signed int j = 0; j < h.nodes.size(); j++) {
				printf("hnodes.id  = %10d\n", h.nodes.at(j).id);
				printf("hnodes.val = %10f\n", h.nodes.at(j).val);
			}
			printf("\n");
			for (signed int j = 0; j < parents.size(); j++) {
				printf("parents.id  = %10d\n", parents.at(j).id);
				printf("parents.val = %10f\n", parents.at(j).val);
			}
		}

		printf("\n========================================================================\n");
		*/
	}

	dc.rowval = subMember(dc.rowval, 0, dc.colptr.at(dc.colptr.size() - 1) - 1);

	for (signed int i = 0; i < dc.rowval.size(); i++) {
		dc.rowval[i].id = dc.revP.at(dc.rowval[i].id);
	}

	result.colptr = dc.colptr;
	result.rowval = dc.rowval;
	result.P = dc.P;
	result.revP = dc.revP;
	result.distances = distances;

	/*
	result.checkInt.resize(parents.size());
	for (signed int i = 0; i < parents.size(); i++) {
		result.checkInt.at(i) = parents.at(i).id;
	}

	result.checkDouble.resize(parents.size());
	for (signed int i = 0; i < parents.size(); i++) {
		result.checkDouble.at(i) = parents.at(i).val;
	}
	*/

	return(result);
}
