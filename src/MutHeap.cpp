// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "rcpplist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

using namespace std;

// struct Node
// {
//   double val; // value
//   signed int id; // id
//   signed int rank; // larger ranks get picked last
//   
//   bool operator>(const Node& other);
//   bool operator>=(const Node& other);
// };

bool Node::operator>(const Node& other) const
{
  if (rank == other.rank) {
    if (val > other.val) {
      return(true);
    }
    else {
      return(false);
    }
  }
  else if (rank > other.rank) {
    return(false);
  }
  else {
    return(true);
  }
  
  /*
   if (val > other.val) { // rank < other.rank &&
   return(true);
   }
   else {
   return(false);
   }
   */
}

bool Node::operator>=(const Node& other) const
{
  if (rank == other.rank) {
    if (val >= other.val) {
      return(true);
    }
    else {
      return(false);
    }
  }
  else if (rank > other.rank) {
    return(false);
  }
  else {
    return(true);
  }
  
  /*
   if ( val >= other.val) { // rank <= other.rank &&
   return(true);
   }
   else {
   return(false);
   }
   */
}

// struct MutHeap
// {
//   vector<Node> nodes; // Vector containing the nodes of the heap
//   vector<signed int> lookup; // Vector providing a lookup for the nodes
// };

void _swap(MutHeap *h, signed int a, signed int b)
{
  h->lookup[h->nodes[a].id] = b;
  h->lookup[h->nodes[b].id] = a;
  
  Node tempNode = h->nodes[a];
  h->nodes[a] = h->nodes[b];
  h->nodes[b] = tempNode;
}

signed int _moveDown(MutHeap *h, signed int hInd)
{
  Node pivot = h->nodes[hInd];
  
  if (2 * hInd + 2 <= h->nodes.size() - 1) { // If both children exist:
    
    if (h->nodes[2 * hInd + 1] >= h->nodes[2 * hInd + 2]) { // If the left child is larger:
      
      if (h->nodes[2 * hInd + 1] >= pivot) { // If the child is larger than the parent:
        _swap(h, hInd, 2 * hInd + 1);
        return(2 * hInd + 1);
      }
      else { // No swap occuring:
        return(h->nodes.size() - 1);
      }
      
    }
    else { // If the right child is larger:
      
      if (h->nodes[2 * hInd + 2] >= pivot) { // If the child is larger than the parent:
        _swap(h, hInd, 2 * hInd + 2);
        return(2 * hInd + 2);
      }
      else { // No swap occuring:
        return(h->nodes.size() - 1);
      }
      
    }
    
  }
  else if (2 * hInd + 1 <= h->nodes.size() - 1) { // If only one child exists:
    
    if (h->nodes[2 * hInd + 1] > pivot) { // If the child is larger than the parent:
      _swap(h, hInd, 2 * hInd + 1);
      return(2 * hInd + 1);
    }
    
  }
  else { // If no children exist:
    
    return(h->nodes.size() - 1);
    
  }
  
  return(h->nodes.size() - 1);
}

Node topNode(MutHeap *h)
{
  return(h->nodes.front());
}

Node topNode_rankUpdate(MutHeap *h)
{
  h->nodes.at(0).rank = numeric_limits<signed int>::max();
  return(h->nodes.front());
}

double update(MutHeap *h, signed int hInd, double hVal)
{
  signed int tempInd = h->lookup[hInd];
  
  if (h->nodes[tempInd].val > hVal) {
    
    h->nodes[tempInd].val = hVal;
    h->nodes[tempInd].id = hInd;
    
    while (tempInd < h->nodes.size() - 1) {
      tempInd = _moveDown(h, tempInd);
    }
    
    return(hVal);
    
  }
  else {
    
    return(h->nodes.at(hInd).val);
    
  }
}