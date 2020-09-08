#pragma once

#include "MutHeap.h"
#include "math.h"
#include "algorithm"
using namespace std;

struct Member
{
  double val;
  signed int id;
  
  bool operator<(const Member& other);
};

bool Member::operator<(const Member& other)
{
  if (val < other.val && id > other.id ) {
    return(true);
  }
  else {
    return(false);
  }
}

bool compareMember(Member const a, Member const b)
{
  if (a.val < b.val) {
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
}

vector<Member> subMember(vector<Member> vec, signed int a, signed int b)
{
  auto first = vec.begin() + a;
  auto last = vec.begin() + b;
  
  vector<Member> subvec(first, last);
  
  return(subvec);
}

void _determineChildren(MutHeap *h, ChildList *dc, vector<Member> *parents, Node pivot, vector<Member> *buffer, double rho, function<double(signed int, signed int)> dist2Func)
{
  double distToParent = parents->at(pivot.id).val;
  double lengthScale = pivot.val;
  signed int iterBuffer = 0;
  
  Member candidate;
  double dist, dist2, newDist;
  vector<Member> viewBuffer;	
  
  signed int start = dc->colptr[dc->revP[parents->at(pivot.id).id]];
  signed int end = dc->colptr[dc->revP[parents->at(pivot.id).id] + 1];
  
  for (signed int i = start; i < end; i++) {
    
    candidate = dc->rowval.at(i);	
    dist2 = dist2Func(candidate.id, pivot.id); 
    
    if (dc->revP.at(candidate.id) == -1 && dist2 <= pow(lengthScale * rho, 2.0)) {
      
      dist = sqrt(dist2);
      
      buffer->at(iterBuffer) = assignMember(dist, candidate.id);
      iterBuffer++;
      
      newDist = update(h, candidate.id, dist);
      
      if (dist + rho * newDist <= rho * lengthScale && dist < parents->at(candidate.id).val) {
        parents->at(candidate.id) = assignMember(dist, pivot.id);
      }
      
    }
    else if (candidate.val > distToParent + lengthScale * rho) {
      break;
    }
  }
  
  viewBuffer = subMember(*buffer, 0, iterBuffer);
  sort(viewBuffer.begin(), viewBuffer.end(), compareMember);
  newParent(dc, pivot.id); // printf("%10d", pivot.id);
  newChildren(dc, viewBuffer);
}	

struct output
{
  vector<signed int> colptr;
  vector<Member> rowval;
  vector<signed int> P;
  vector<signed int> revP;
  vector<double> distances;
};

output sortSparse(signed int N, double rho, function<double(signed int, signed int)> dist2Func, signed int initInd)
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
  }
  
  dc.rowval = subMember(dc.rowval, 0, dc.colptr.at(dc.colptr.size() - 1));
  
  for (signed int i = 0; i < dc.rowval.size(); i++) {
    dc.rowval[i].id = dc.revP.at(dc.rowval[i].id);
  }
  
  result.colptr = dc.colptr;
  result.rowval = dc.rowval;
  result.P = dc.P;
  result.revP = dc.revP;
  result.distances = distances;
  
  return(result);
}

output predSortSparse(signed int NTrain, signed int NTest, double rho, function<double(signed int, signed int)> dist2Func, signed int initInd)
{
  signed int N = NTrain + NTest;
  
  MutHeap h;
  ChildList dc;
  vector<Member> nodeBuffer;
  vector<double> distances;
  
  vector<Member> viewBuffer;
  vector<Member> parents;
  
  output result;
  
  h.nodes.resize(N);
  h.lookup.resize(N);
  
  for (signed int i = 0; i < NTrain; i++) {
    h.nodes[i].val = numeric_limits<double>::max();
    h.nodes[i].id = i;
    h.lookup[i] = i;
    
    h.nodes[i].rank = 0;
  }
  
  for (signed int i = NTrain; i < N; i++) {
    h.nodes[i].val = numeric_limits<double>::max();
    h.nodes[i].id = i;
    h.lookup[i] = i;
    
    h.nodes[i].rank = 1;
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
  }
  
  dc.rowval = subMember(dc.rowval, 0, dc.colptr.at(dc.colptr.size() - 1));
  
  for (signed int i = 0; i < dc.rowval.size(); i++) {
    dc.rowval[i].id = dc.revP.at(dc.rowval[i].id);
  }
  
  result.colptr = dc.colptr;
  result.rowval = dc.rowval;
  result.P = dc.P;
  result.revP = dc.revP;
  result.distances = distances;
  
  return(result);
}
