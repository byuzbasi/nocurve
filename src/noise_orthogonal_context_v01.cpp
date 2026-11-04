// Isolated implementation of the frozen NOC v01 specification.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List noc_basis_cpp(const arma::mat& B, const arma::vec& a, bool project) {
  arma::mat R = B;
  if (project && B.n_rows) {
    double aa = arma::dot(a,a);
    if (!(aa > 0)) stop("Nonpositive target norm");
    R -= (B*a)*a.t()/aa;
  }
  arma::mat U,V; arma::vec s;
  if (!R.n_rows) return List::create(_["Q"]=arma::mat(a.n_elem,0),
    _["singular"]=s, _["rank"]=0, _["cutoff"]=1e-12);
  if (!arma::svd_econ(U,s,V,R.t(),"both","std")) stop("Context SVD failed");
  double cut = 1e-12*std::max(1.0,s.max());
  if (arma::any((s > cut/10) % (s < cut*10))) stop("Ambiguous context rank");
  arma::uvec keep=arma::find(s>cut);
  arma::mat Q=U.cols(keep);
  return List::create(_["Q"]=Q, _["singular"]=s,
    _["rank"]=keep.n_elem, _["cutoff"]=cut);
}

static double lse(const arma::vec& x) {
  double m=x.max();
  return m+std::log(arma::sum(arma::exp(x-m)));
}

// Stable Gaussian-mixture posterior and the two marginal log densities.
// [[Rcpp::export]]
List noc_posterior_cpp(const arma::vec& d, const arma::vec& v,
                      const arma::vec& pi, const arma::vec& sd,
                      const arma::vec& c, double omega) {
  if (d.n_elem!=v.n_elem || d.n_elem!=c.n_elem || pi.n_elem!=sd.n_elem ||
      omega<0 || omega>1 || !arma::all(v>0) || !arma::all(pi>=0) ||
      !arma::all(sd>=0) || std::abs(arma::sum(pi)-1)>1e-8)
    stop("Invalid mixture input");
  int n=d.n_elem, K=pi.n_elem;
  double p0=arma::sum(pi.elem(arma::find(sd==0)));
  arma::vec mean(n),factor(n),log0(n),log1(n);
  arma::mat prob(n,K),weights(n,K);
  for (int i=0;i<n;++i) {
    if (!std::isfinite(c[i]) || c[i]<-1 || c[i]>1) stop("Invalid context score");
    double tilt=std::pow(4.0,c[i]), Z=p0+tilt*(1-p0);
    arma::vec lp(K), l0(K), l1(K), sh(K), pw(K);
    for (int k=0;k<K;++k) {
      double p1=pi[k]*(sd[k]==0?1:tilt)/Z;
      pw[k]=(1-omega)*pi[k]+omega*p1;
      double vv=v[i]+sd[k]*sd[k];
      double lognormal=-0.5*(std::log(2*arma::datum::pi*vv)+d[i]*d[i]/vv);
      l0[k]=std::log(pi[k])+lognormal;
      l1[k]=std::log(p1)+lognormal;
      lp[k]=std::log(pw[k])+lognormal;
      sh[k]=sd[k]*sd[k]/vv;
    }
    log0[i]=lse(l0); log1[i]=lse(l1);
    arma::vec post=arma::exp(lp-lse(lp));
    factor[i]=arma::dot(post,sh); mean[i]=d[i]*factor[i];
    prob.row(i)=post.t(); weights.row(i)=pw.t();
  }
  return List::create(_["mean"]=mean, _["factor"]=factor,
    _["log0"]=log0, _["log1"]=log1,
    _["posterior_prob"]=prob, _["prior_weights"]=weights);
}
