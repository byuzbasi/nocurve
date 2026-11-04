# Selected implementation from R/noc_compact_v01.R.
ncp_yhash<-function(y)digest::digest(y,algo='sha256')

ncp_impulse<-function(n) {
  N<-2L*n;J<-as.integer(log2(N));stopifnot(2^J==N)
  e<-numeric(N);e[1]<-1;w<-noc_wd(e)
  g<-matrix(noc_details(w),N,J)
  for(j in 0:(J-1L))w<-wavethresh::putD(w,j,numeric(N))
  s<-as.numeric(wavethresh::AvBasis(wavethresh::convert(w)))
  list(n=n,N=N,J=J,g=g,scaling_impulse=s,format='NOC-level-stream-v01')
}

ncp_descriptor<-function(y,I,levels) {
  stopifnot(length(levels)==I$J)
  for(j in 0:(I$J-1L))stopifnot(levels[[j+1L]]$valid,levels[[j+1L]]$j==j,
    levels[[j+1L]]$n==I$n,levels[[j+1L]]$response_sha256==ncp_yhash(y))
  bind<-function(name)unlist(lapply(levels,`[[`,name),use.names=FALSE)
  z<-c(I,list(aa=bind('aa'),native=bind('native'),active=bind('active'),
    analysis=bind('analysis'),parent=bind('parent'),q=bind('q'),qraw=bind('qraw'),
    neighbors=do.call(rbind,lapply(levels,`[[`,'neighbors')),
    energy=bind('energy'),raw_energy=bind('raw_energy'),response_sha256=ncp_yhash(y)))
  stopifnot(length(z$aa)==I$J*I$N,identical(as.integer(z$q),z$q),
    identical(as.integer(z$qraw),z$qraw))
  z
}

ncp_analysis<-function(D,y) {
  stopifnot(length(y)==D$n,ncp_yhash(y)==D$response_sha256)
  D$analysis
}

ncp_energy<-function(D,y,raw) {
  stopifnot(ncp_yhash(y)==D$response_sha256)
  if(raw)D$raw_energy else D$energy
}

ncp_rewrite<-function(fun,replacements,env) {
  counts<-integer(length(replacements))
  walk<-function(e) {
    for(i in seq_along(replacements))if(identical(e,replacements[[i]][[1]])) {
      counts[i]<<-counts[i]+1L;return(replacements[[i]][[2]])
    }
    if(is.call(e))for(k in seq_along(e))e[[k]]<-walk(e[[k]])
    e
  }
  body(fun)<-walk(body(fun));stopifnot(all(counts==1L));environment(fun)<-env;fun
}

ncp_functions<-function() {
  env<-new.env(parent=environment(ncp_functions))
  env$ncp_synthesis<-n3a_synthesis
  env$noc_context<-ncp_rewrite(noc_context,list(list(
    quote(vapply(Q,function(u)sum(crossprod(u,y)^2),0.0)),quote(ncp_energy(D,y,raw)))),env)
  env$noc_base<-ncp_rewrite(noc_base,list(list(quote(drop(D$A%*%y)),
    quote(ncp_analysis(D,y)))),env)
  # The original posterior and both row-sum gates still run at every level.
  # Retain a certificate instead of all two dense probability matrices at once.
  env$noc_evidence<-ncp_rewrite(noc_evidence,list(list(quote(levels[[j]]<-z),
    quote(levels[[j]]<-list(storage='validated-level-summary',
      observations=length(ix),components=length(p$pi))))),env)
  env$noc_output<-noc_output;environment(env$noc_output)<-env
  env$noc_validate_output<-ncp_rewrite(noc_validate_output,list(list(
    quote(drop(D$synthesis%*%e$theta+D$scaling%*%base$y)),
    quote(ncp_synthesis(D,e$theta,base$y)))),env)
  env$niv_check_base<-ncp_rewrite(niv_check_base,list(list(quote(drop(D$A%*%base$y)),
    quote(ncp_analysis(D,base$y)))),env)
  env
}

# Selected implementation from R/noc_three_real_adapters_v01.R.
n3a_extension <- function(n) {
  stopifnot(length(n)==1L,is.finite(n),n==as.integer(n),n>=4L)
  # Same source indices as mc_extension, omitting only its dense one-hot matrix.
  n<-as.integer(n);m<-as.integer(2^ceiling(log2(n)))
  pad<-c(seq_len(n),if(m>n)seq.int(n,2L*n-m+1L) else integer())
  list(n=n,m=m,N=2L*m,pad_index=pad,index=c(pad,rev(pad)))
}

n3a_impulse <- function(n) {
  ex<-n3a_extension(n);I<-ncp_impulse(ex$m)
  I$n<-ex$n;I$m<-ex$m;I$index<-ex$index;I$pad_index<-ex$pad_index
  I$format<-'NOC-source-copy-stream-v01';I
}

n3a_rows <- function(I,j) {
  a<-ncp_analysis_level_cpp(I$g[,j+1L],I$m)
  if(I$n!=I$m) {
    # Folding columns composes the original linear analysis with the source-copy map.
    a$A<-t(rowsum(t(a$A),I$pad_index,reorder=FALSE));dimnames(a$A)<-NULL
    a$aa<-rowSums(a$A^2)
    if(any(sqrt(a$aa[a$aa>0])<=1e-10))stop('Ambiguous near-zero analysis row')
  }
  a
}

n3a_level <- function(y,I,j) {
  a<-n3a_rows(I,j)
  if(j>0L) {
    prev<-n3a_rows(I,j-1L);parent<-as.integer(ncp_parents_cpp(a$A,prev$A,a$aa,prev$aa))
  } else {prev<-list(A=matrix(0,0,0),aa=numeric());parent<-rep(-1L,I$N)}
  z<-ncp_context_level_cpp(a$A,prev$A,a$aa,prev$aa,parent,y,2^(I$J-1L-j))
  dy<-drop(a$A%*%y)
  stopifnot(all(is.finite(c(dy,a$aa,a$native,z$energy,z$raw_energy))),
    all(z$energy>=0),all(z$raw_energy>=0),z$max_orthogonality<=1e-8)
  list(n=I$n,j=j,aa=as.numeric(a$aa),native=as.numeric(a$native),active=a$aa>0,
    analysis=dy,parent=ifelse(parent<0L,NA_integer_,(j-1L)*I$N+parent+1L),
    neighbors=z$neighbors+j*I$N,q=as.integer(z$q),qraw=as.integer(z$qraw),
    energy=as.numeric(z$energy),raw_energy=as.numeric(z$raw_energy),
    singular=z$singular,raw_singular=z$raw_singular,max_orthogonality=z$max_orthogonality,
    response_sha256=ncp_yhash(y),valid=TRUE)
}

n3a_synthesis <- function(D,theta,y) {
  stopifnot(length(theta)==D$N*D$J,length(y)==D$n,all(is.finite(c(theta,y))))
  p<-fft(fft(D$scaling_impulse)*fft(y[D$index]),inverse=TRUE)/D$N
  for(j in 0:(D$J-1L))p<-p+2^(j-D$J)*fft(Conj(fft(D$g[,j+1L]))*
    fft(theta[j*D$N+seq_len(D$N)]),inverse=TRUE)/D$N
  stopifnot(max(abs(Im(p)))<=1e-8*max(1,abs(p)));Re(p)[seq_len(D$n)]
}

