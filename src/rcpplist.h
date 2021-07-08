#ifndef RCPPLIST_H
#define RCPPLIST_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;

// example.cpp

arma::vec fun_Rcpp_example(const double & a, const arma::vec & x);

// covftns.cpp

double cov_expo_iso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

double cov_expo_aniso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

double cov_expo_spacetime_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

double cov_expo_squared_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

double cov_matern_iso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

double cov_matern_aniso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

double cov_matern_scaledim_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

double cov_matern_spacetime_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);

arma::mat cov_matern_ns_forR(const arma::mat & locs1, const arma::mat & locs2, Rcpp::Function sigma, Rcpp::Function smoothness, Rcpp::Function kernel);

arma::mat cal_matern_ns(const arma::mat & locs, const arma::rowvec & sigvec, const arma::rowvec & nuvec, const Rcpp::List& kerList);

double cov_matern_ns_cpp(const arma::rowvec & x1, const arma::rowvec & x2, Rcpp::Function sigma, Rcpp::Function smoothness, Rcpp::Function kernel);

double cov_latentDim_biv_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec covparms);

double cov_latentDim_triv_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec covparms);

double GpGp_matern_spacetime_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec covparms);

// MutHeap.cpp

struct Node
{
  double val; // value
  signed int id; // id
  signed int rank; // larger ranks get picked last

  bool operator>(const Node& other) const;
  bool operator>=(const Node& other) const;
};

// bool Node::operator>(const Node& other);

// bool Node::operator>=(const Node& other);

struct MutHeap
{
  vector<Node> nodes; // Vector containing the nodes of the heap
  vector<signed int> lookup; // Vector providing a lookup for the nodes
};

void _swap(MutHeap *h, signed int a, signed int b);

signed int _moveDown(MutHeap *h, signed int hInd);

Node topNode(MutHeap *h);

Node topNode_rankUpdate(MutHeap *h);

double update(MutHeap *h, signed int hInd, double hVal);

// SortSparse.cpp

struct Member
{
  double val;
  signed int id;

  bool operator<(const Member& other) const;
};

// bool Member::operator<(const Member& other);

bool compareMember(Member const a, Member const b);

Member assignMember(double val, signed int id);

struct ChildList
{
  signed int NParents = 0, NChildren = 0, NBuffer = 0;

  vector<signed int> P; // This array gives for contains the ordering.The i-th parent in the daycare has id P[i]
  vector<signed int> revP; // This array contains as the i-th element the number that the ith parent has with respect to the multiresolution ordering.
  vector<signed int> colptr; // The array that contains the first "child" for every parent
  vector<Member> rowval; // The array that contains the global id-s of the children
};

void newParent(ChildList *dc, signed int idParent);

void newChildren(ChildList *dc, vector<Member> children);

vector<Member> subMember(vector<Member> vec, signed int a, signed int b);

void _determineChildren(MutHeap *h, ChildList *dc, vector<Member> *parents, Node pivot, vector<Member> *buffer, double rho, function<double(signed int, signed int)> dist2Func);

struct output
{
  vector<signed int> colptr;
  vector<Member> rowval;
  vector<signed int> P;
  vector<signed int> revP;
  vector<double> distances;
};

output sortSparse(signed int N, double rho, function<double(signed int, signed int)> dist2Func, signed int initInd);

output predSortSparse(signed int NTrain, signed int NTest, double rho, function<double(signed int, signed int)> dist2Func, signed int initInd);

// ic0.cpp

Rcpp::NumericVector mat2vals(const Rcpp::NumericVector ptrs, const Rcpp::NumericVector inds, const Rcpp::NumericMatrix covmat);

#endif
