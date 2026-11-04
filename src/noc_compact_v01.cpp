// [[Rcpp::depends(RcppArmadillo)]]
// Reuse the original SVD and posterior implementation verbatim, not a Gram surrogate.
#include "noise_orthogonal_context_v01.h"
using Rcpp::List;

// [[Rcpp::export]]
Rcpp::List ncp_analysis_level_cpp(const arma::vec& g, int n) {
  int N=2*n;
  if (g.n_elem!=static_cast<unsigned>(N)) Rcpp::stop("Impulse length mismatch");
  arma::mat A(N,n);
  for(int l=0;l<n;++l) for(int k=0;k<N;++k)
    A(k,l)=g[(k-l+N)%N]+g[(k+1+l)%N];
  arma::vec aa(N);
  // R rowSums accumulates double products in extended precision.
  for(int k=0;k<N;++k) {
    long double sum=0;for(int l=0;l<n;++l){double v=A(k,l)*A(k,l);sum+=v;}
    aa[k]=static_cast<double>(sum);
  }
  double native=arma::dot(g,g);
  for(int i=0;i<N;++i) if(aa[i]>0 && std::sqrt(aa[i])<=1e-10)
    Rcpp::stop("Ambiguous near-zero analysis row");
  return Rcpp::List::create(Rcpp::Named("A")=A,Rcpp::Named("aa")=aa,
    Rcpp::Named("native")=arma::vec(N,arma::fill::value(native)));
}

// [[Rcpp::export]]
arma::ivec ncp_parents_cpp(const arma::mat& A,const arma::mat& previous,
  const arma::vec& aa,const arma::vec& previous_aa,int block=64) {
  int N=A.n_rows;
  if(previous.n_rows!=A.n_rows || previous.n_cols!=A.n_cols || block<1)
    Rcpp::stop("Parent geometry dimensions differ");
  arma::ivec parent(N);parent.fill(-1);
  arma::mat ps=arma::square(previous);
  for(int first=0;first<N;first+=block) {
    Rcpp::checkUserInterrupt();int last=std::min(first+block-1,N-1);
    arma::mat overlap=arma::square(A.rows(first,last))*ps.t();
    for(int i=first;i<=last;++i) if(aa[i]>0) {
      double best=-arma::datum::inf;
      for(int p=0;p<N;++p) if(previous_aa[p]>0) {
        overlap(i-first,p)/=aa[i]*previous_aa[p];
        best=std::max(best,overlap(i-first,p));
      }
      // Same absolute tie rule and recorded index order as noc_design.
      for(int p=0;p<N;++p) if(previous_aa[p]>0 &&
        overlap(i-first,p)>=best-1e-12) {parent[i]=p;break;}
    }
  }
  return parent;
}

// [[Rcpp::export]]
Rcpp::List ncp_context_level_cpp(const arma::mat& A,const arma::mat& previous,
  const arma::vec& aa,const arma::vec& previous_aa,const arma::ivec& parent,
  const arma::vec& y,int h) {
  int N=A.n_rows,n=A.n_cols;
  if(y.n_elem!=static_cast<unsigned>(n) || parent.n_elem!=static_cast<unsigned>(N))
    Rcpp::stop("Context dimensions differ");
  arma::vec energy(N,arma::fill::zeros),raw_energy(N,arma::fill::zeros);
  arma::ivec rank(N,arma::fill::zeros),raw_rank(N,arma::fill::zeros);
  arma::mat singular(N,3),raw_singular(N,3);singular.fill(NA_REAL);raw_singular.fill(NA_REAL);
  arma::imat near(N,2);double orth=0;
  for(int k=0;k<N;++k) {
    if(k%256==0)Rcpp::checkUserInterrupt();
    int ids[2]={(k-h+N)%N,(k+h)%N};near(k,0)=ids[0]+1;near(k,1)=ids[1]+1;
    int rows=0;for(int j=0;j<2;++j)if(aa[ids[j]]>0)++rows;
    if(parent[k]>=0)++rows;
    if(aa[k]<=0)continue;
    arma::mat B(rows,n);int r=0;
    for(int j=0;j<2;++j)if(aa[ids[j]]>0)B.row(r++)=A.row(ids[j])/std::sqrt(aa[ids[j]]);
    if(parent[k]>=0)B.row(r++)=previous.row(parent[k])/std::sqrt(previous_aa[parent[k]]);
    arma::vec a=A.row(k).t();
    Rcpp::List projected=noc_basis_cpp(B,a,true),raw=noc_basis_cpp(B,a,false);
    arma::mat Q=Rcpp::as<arma::mat>(projected["Q"]),Qr=Rcpp::as<arma::mat>(raw["Q"]);
    arma::vec s=Rcpp::as<arma::vec>(projected["singular"]),sr=Rcpp::as<arma::vec>(raw["singular"]);
    rank[k]=Rcpp::as<int>(projected["rank"]);raw_rank[k]=Rcpp::as<int>(raw["rank"]);
    energy[k]=arma::dot(Q.t()*y,Q.t()*y);raw_energy[k]=arma::dot(Qr.t()*y,Qr.t()*y);
    for(unsigned j=0;j<s.n_elem;++j)singular(k,j)=s[j];
    for(unsigned j=0;j<sr.n_elem;++j)raw_singular(k,j)=sr[j];
    if(Q.n_cols)orth=std::max(orth,arma::abs(a.t()*Q).max()/std::sqrt(aa[k]));
  }
  return Rcpp::List::create(Rcpp::Named("energy")=energy,Rcpp::Named("raw_energy")=raw_energy,
    Rcpp::Named("q")=rank,Rcpp::Named("qraw")=raw_rank,Rcpp::Named("neighbors")=near,
    Rcpp::Named("singular")=singular,Rcpp::Named("raw_singular")=raw_singular,
    Rcpp::Named("max_orthogonality")=orth);
}
