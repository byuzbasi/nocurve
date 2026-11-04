# Selected implementation from R/noise_orthogonal_context_v01.R.
noc_methods <- c("EB0","NOC","RAW")

noc_ash_options <- list(pointmass=TRUE,prior="nullbiased",nullweight=10,
  gridmult=64,mixcompdist="normal",outputlevel=2,fixg=FALSE)

noc_assert <- function(a,b) {
  stopifnot(length(a)==length(b),all(is.finite(a)),all(is.finite(b)),
    max(c(0,abs(a-b)))<=1e-8*max(1,abs(a),abs(b)))
  invisible(TRUE)
}

noc_details <- function(w) unlist(lapply(0:(w$nlevels-1L),
  function(j) wavethresh::accessD(w,j)),use.names=FALSE)

noc_wd <- function(y) wavethresh::wd(y,type="station",family="DaubExPhase",filter.number=4)

noc_inverse <- function(w,theta,n) {
  N<-2L*n
  for(j in 0:(w$nlevels-1L))w<-wavethresh::putD(w,j,theta[j*N+seq_len(N)])
  as.numeric(wavethresh::AvBasis(wavethresh::convert(w)))[seq_len(n)]
}

noc_design <- function(n=128L) {
  N<-2L*n;J<-as.integer(log2(N));stopifnot(2^J==N)
  E<-rbind(diag(n),diag(n)[n:1,])
  W<-vapply(seq_len(N),function(k)noc_details(noc_wd(diag(N)[,k])),numeric(J*N))
  A<-W%*%E;aa<-rowSums(A^2);native<-rowSums(W^2)
  active<-aa>0
  if(any(sqrt(aa[active])<=1e-10))stop("Ambiguous near-zero analysis row")
  Q<-Qraw<-Blist<-singular<-vector("list",J*N)
  q<-qr<-integer(J*N);parent<-rep(NA_integer_,J*N)
  neighbors<-matrix(NA_integer_,J*N,2)
  for(j in 0:(J-1L)) {
    ix<-j*N+seq_len(N);h<-2^(J-1L-j)
    if(j>0) {
      prev<-(j-1L)*N+seq_len(N);valid<-which(active[prev])
      overlap<-(A[ix,,drop=FALSE]^2)%*%t(A[prev,,drop=FALSE]^2)
      overlap<-overlap/outer(aa[ix],aa[prev])
    }
    for(k in seq_len(N)) {
      i<-ix[k];near<-j*N+((c(k-h,k+h)-1L)%%N+1L);neighbors[i,]<-near
      if(j>0 && active[i] && length(valid)) {
        vals<-overlap[k,valid];pick<-valid[which(vals>=max(vals)-1e-12)[1]]
        parent[i]<-prev[pick];near<-c(near,parent[i])
      }
      near<-near[active[near]]
      B<-A[near,,drop=FALSE]/sqrt(aa[near]);Blist[[i]]<-B
      if(!active[i]) {Q[[i]]<-Qraw[[i]]<-matrix(0,n,0);next}
      z<-noc_basis_cpp(B,A[i,],TRUE);raw<-noc_basis_cpp(B,A[i,],FALSE)
      Q[[i]]<-z$Q;Qraw[[i]]<-raw$Q;q[i]<-z$rank;qr[i]<-raw$rank
      singular[[i]]<-list(projected=z$singular,raw=raw$singular)
    }
  }
  zero<-noc_wd(rep(0,N));H<-J*N
  synthesis<-vapply(seq_len(H),function(i){v<-numeric(H);v[i]<-1;noc_inverse(zero,v,n)},numeric(n))
  scaling<-vapply(seq_len(n),function(i)noc_inverse(noc_wd(E[,i]),numeric(H),n),numeric(n))
  list(n=n,N=N,J=J,E=E,W=W,A=A,aa=aa,native=native,active=active,
    synthesis=synthesis,scaling=scaling,
    Q=Q,Qraw=Qraw,B=Blist,q=q,qraw=qr,parent=parent,neighbors=neighbors,singular=singular)
}

noc_context <- function(y,sigma,D,raw=FALSE) {
  Q<-if(raw)D$Qraw else D$Q;q<-if(raw)D$qraw else D$q
  energy<-vapply(Q,function(u)sum(crossprod(u,y)^2),0.0)
  c<-rep(0,length(q));ix<-q>0
  c[ix]<-2*pchisq(energy[ix]/sigma^2,df=q[ix])-1
  list(c=c,energy=energy,q=q,raw=raw)
}

noc_ash_fit <- function(d,se) {
  .GlobalEnv$.noc_native<-list()
  suppressMessages(trace("estimate_mixprop",where=asNamespace("ashr"),print=FALSE,
    exit=quote({.a<-returnValue();.GlobalEnv$.noc_native[[length(.GlobalEnv$.noc_native)+1L]]<-
      list(optreturn=.a$optreturn,optmethod=.a$optmethod,control=control)})))
  on.exit(suppressMessages(untrace("estimate_mixprop",where=asNamespace("ashr"))))
  warnings<-character()
  f<-withCallingHandlers(do.call(ashr::ash,c(list(betahat=d,sebetahat=se),noc_ash_options)),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  native<-.GlobalEnv$.noc_native;g<-f$fitted_g
  stopifnot(!length(warnings),length(native)>0L,
    all(vapply(native,function(z)isTRUE(z$optreturn$converged),TRUE)),
    inherits(g,"normalmix"),all(g$mean==0),all(is.finite(c(g$pi,g$sd))),
    all(g$pi>=0),all(g$sd>=0),abs(sum(g$pi)-1)<=1e-8)
  prior<-list(pi=c(sum(g$pi[g$sd==0]),g$pi[g$sd>0]),sd=c(0,g$sd[g$sd>0]))
  list(fit=f,prior=prior,native=native,warnings=warnings,options=noc_ash_options)
}

noc_base <- function(y,sigma,D,native_variance=FALSE) {
  stopifnot(length(y)==D$n,all(is.finite(y)),length(sigma)==1,is.finite(sigma),sigma>=0)
  if(sigma==0) {
    if(!all(y==y[1]))stop("Zero estimated noise on nonconstant data is unsupported")
    return(list(constant=y[1],y=y,sigma=sigma))
  }
  w<-noc_wd(c(y,rev(y)));d<-noc_details(w);noc_assert(d,drop(D$A%*%y))
  v<-sigma^2*if(native_variance)D$native else D$aa
  fits<-vector("list",D$J);pm<-numeric(length(d));began<-proc.time()[[3]]
  for(j in 0:(D$J-1L)) {
    ix<-j*D$N+seq_len(D$N);ix<-ix[D$active[ix]]
    if(!length(ix))next
    f<-noc_ash_fit(d[ix],sqrt(v[ix]));f$index<-ix;fits[[j+1L]]<-f
    z<-noc_posterior_cpp(d[ix],v[ix],f$prior$pi,f$prior$sd,rep(0,length(ix)),0)
    noc_assert(z$mean,ashr::get_pm(f$fit));pm[ix]<-z$mean
  }
  list(y=y,sigma=sigma,w=w,d=d,v=v,fits=fits,pm=pm,native_variance=native_variance,
    seconds=proc.time()[[3]]-began)
}

noc_score <- function(omega,log0,log1,n) {
  H<-length(log0);if(!H)return(c(value=-omega^2/2,derivative=-omega,curvature=-1,kkt=abs(omega)))
  r<-expm1(log1-log0);z<-1+omega*r;t<-r/z
  derivative<-n/H*sum(t)-omega
  # Endpoint KKT uses the appropriate feasible direction.
  violation<-if(omega==0)max(0,derivative) else if(omega==1)max(0,-derivative) else abs(derivative)
  c(value=n/H*sum(log0+log1p(omega*r))-omega^2/2,
    derivative=derivative,curvature=-n/H*sum(t^2)-1,
    kkt=violation/max(1,n/H*sum(abs(t))+abs(omega)))
}

noc_optimize <- function(log0,log1,n) {
  lo<-0;hi<-1;iterations<-0L
  if(noc_score(0,log0,log1,n)["derivative"]<=0)omega<-0
  else if(noc_score(1,log0,log1,n)["derivative"]>=0)omega<-1
  else {
    repeat {
      omega<-(lo+hi)/2;iterations<-iterations+1L
      z<-noc_score(omega,log0,log1,n)
      if(z["derivative"]>0)lo<-omega else hi<-omega
      if(hi-lo<=1e-10 && z["kkt"]<=1e-8)break
      if(iterations>=80L)stop("Scalar optimizer did not meet frozen gates")
    }
  }
  z<-noc_score(omega,log0,log1,n);stopifnot(z["kkt"]<=1e-8)
  list(omega=omega,score=z,bracket=c(lo,hi),iterations=iterations,
    endpoints=rbind(noc_score(0,log0,log1,n),noc_score(1,log0,log1,n)))
}

noc_evidence <- function(base,c,D,omega) {
  theta<-factor<-numeric(length(base$d));log0<-log1<-rep(NA_real_,length(base$d));levels<-vector("list",D$J)
  for(j in seq_len(D$J)) {
    f<-base$fits[[j]];if(is.null(f))next
    ix<-f$index;p<-f$prior
    z<-noc_posterior_cpp(base$d[ix],base$v[ix],p$pi,p$sd,c[ix],omega)
    theta[ix]<-z$mean;factor[ix]<-z$factor;log0[ix]<-z$log0;log1[ix]<-z$log1;levels[[j]]<-z
    stopifnot(all(z$factor>=0),all(z$factor<=1+1e-12),all(is.finite(z$mean)))
    noc_assert(rowSums(z$posterior_prob),rep(1,length(ix)))
    noc_assert(rowSums(z$prior_weights),rep(1,length(ix)))
  }
  list(theta=theta,factor=factor,log0=log0,log1=log1,levels=levels)
}

noc_output <- function(base,D,method) {
  stopifnot(method%in%noc_methods);started<-proc.time()[[3]]
  if(!is.null(base$constant))return(list(prediction=rep(base$constant,D$n),method=method,
    constant=TRUE,omega=0,valid=TRUE,seconds=proc.time()[[3]]-started))
  context<-if(method=="EB0")list(c=rep(0,length(D$aa)),q=rep(0,length(D$aa))) else
    noc_context(base$y,base$sigma,D,raw=method=="RAW")
  z<-noc_evidence(base,context$c,D,0)
  opt<-if(method=="EB0")list(omega=0,score=noc_score(0,z$log0[D$active],z$log1[D$active],D$n),
    iterations=0L,bracket=c(0,0)) else noc_optimize(z$log0[D$active],z$log1[D$active],D$n)
  if(opt$omega!=0)z<-noc_evidence(base,context$c,D,opt$omega)
  if(opt$omega==0)noc_assert(z$theta,base$pm)
  p<-noc_inverse(base$w,z$theta,D$n) # Exactly one new complete curve per call.
  stopifnot(length(p)==D$n,all(is.finite(p)),all(z$theta*base$d>=-1e-14))
  list(prediction=p,method=method,context=context,omega=opt$omega,optimization=opt,
    posterior=z,valid=TRUE,seconds=proc.time()[[3]]-started)
}

noc_validate_output <- function(z,base,D,method) {
  stopifnot(identical(z$method,method),isTRUE(z$valid),length(z$prediction)==D$n,
    all(is.finite(z$prediction)),z$omega>=0,z$omega<=1)
  if(isTRUE(z$constant)){noc_assert(z$prediction,base$y);return(invisible(TRUE))}
  context<-if(method=="EB0")list(c=rep(0,length(D$aa))) else noc_context(base$y,base$sigma,D,method=="RAW")
  noc_assert(z$context$c,context$c)
  e<-noc_evidence(base,context$c,D,z$omega)
  noc_assert(z$posterior$theta,e$theta);noc_assert(z$posterior$log0[D$active],e$log0[D$active])
  noc_assert(z$posterior$log1[D$active],e$log1[D$active])
  sc<-noc_score(z$omega,e$log0[D$active],e$log1[D$active],D$n)
  stopifnot(sc["kkt"]<=1e-8)
  # Check the same saved output by an independent linear synthesis matrix;
  # this is validation, not another parameter setting or fitted curve.
  noc_assert(z$prediction,drop(D$synthesis%*%e$theta+D$scaling%*%base$y))
  invisible(TRUE)
}

# Selected implementation from R/multiscale_aggregation_v01.R.
ms_sigma <- function(y) {d<-diff(y);median(abs(d-median(d)))/(sqrt(2)*qnorm(.75))}

# Selected implementation from R/noc_independent_validation_v01.R.
niv_pm <- function(d,v,pi,sd) {
  stopifnot(length(d)==length(v),all(v>0),all(pi>=0),abs(sum(pi)-1)<=1e-8)
  lp<-matrix(vapply(seq_along(pi),function(k)log(pi[k])+dnorm(d,sd=sqrt(v+sd[k]^2),log=TRUE),numeric(length(d))),nrow=length(d))
  top<-apply(lp,1,max);pr<-exp(lp-top);pr<-pr/rowSums(pr)
  sh<-matrix(vapply(sd,function(s)s^2/(v+s^2),numeric(length(d))),nrow=length(d))
  d*rowSums(pr*sh)
}

niv_check_native <- function(f) {
  stopifnot(length(f)>0L,all(vapply(f,function(a)isTRUE(a$optreturn$converged),TRUE)))
}

niv_check_base <- function(base,D) {
  stopifnot(!base$native_variance,length(base$y)==D$n,base$sigma>0)
  noc_assert(base$d,drop(D$A%*%base$y));noc_assert(base$v,base$sigma^2*D$aa)
  for(f in base$fits)if(!is.null(f)) {
    niv_check_native(f$native);stopifnot(!length(f$warnings),identical(f$options,noc_ash_options))
    p<-f$prior;stopifnot(all(p$sd>=0),all(p$pi>=0),abs(sum(p$pi)-1)<=1e-8)
    noc_assert(base$pm[f$index],niv_pm(base$d[f$index],base$v[f$index],p$pi,p$sd))
  }
  invisible(TRUE)
}

