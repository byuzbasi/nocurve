nt_b<-function(u){z<-numeric(length(u));k<-abs(u)<1;z[k]<-(1-u[k]^2)^3;z}

nt_w<-function(u){z<-numeric(length(u));k<-u>0&u<1;z[k]<-sin(pi*u[k])^4;z}

nt_s<-function(u){z<-numeric(length(u));z[u>=1]<-1;k<-u>0&u<1
  z[k]<-10*u[k]^3-15*u[k]^4+6*u[k]^5;z}

nt_eval_mean<-function(expression,x) {
  stopifnot(is.character(expression),length(expression)==1L,all(is.finite(x)))
  e<-parse(text=expression);stopifnot(length(e)==1L)
  allowed<-c('x','pi','sin','cos','b','w','s','+','-','*','/','^','(')
  stopifnot(all(all.names(e)%in%allowed))
  z<-eval(e[[1]],envir=list(x=x,b=nt_b,w=nt_w,s=nt_s),enclos=baseenv())
  stopifnot(length(z)==length(x),all(is.finite(z)));as.numeric(z)
}

niv_options <- function() {
  p<-getFromNamespace('setAshParam.gaus','smashr')(list())
  stopifnot(p$gridmult==2) # Public helper default; native mean fitting overrides it.
  p$gridmult<-64
  stopifnot(p$gridmult==64,p$pointmass,p$prior=='nullbiased',p$nullweight==10,
            p$mixcompdist=='normal',p$outputlevel==2,!p$fixg)
  p
}

niv_smash <- function(y,sigma) {
  expected<-niv_options()
  .GlobalEnv$.niv_native<-list();.GlobalEnv$.niv_smash_records<-list()
  suppressMessages(trace('estimate_mixprop',where=asNamespace('ashr'),print=FALSE,
    exit=quote({.v<-returnValue();.GlobalEnv$.niv_native[[length(.GlobalEnv$.niv_native)+1L]]<-
      list(optreturn=.v$optreturn,optmethod=.v$optmethod,control=control)})))
  on.exit(suppressMessages(untrace('estimate_mixprop',where=asNamespace('ashr'))),add=TRUE)
  suppressMessages(trace('shrink.wc',where=asNamespace('smashr'),print=FALSE,
    tracer=quote({.niv_begin<-length(.GlobalEnv$.niv_native)}),
    exit=quote({.f<-returnValue();.new<-length(.GlobalEnv$.niv_native)-.niv_begin
      .GlobalEnv$.niv_smash_records[[length(.GlobalEnv$.niv_smash_records)+1L]]<-
        list(wc=wc,se=wc.var.sqrt,options=ashparam,fit=.f,
             native=if(.new>0)tail(.GlobalEnv$.niv_native,.new) else list())})))
  on.exit(suppressMessages(untrace('shrink.wc',where=asNamespace('smashr'))),add=TRUE)
  warnings<-character();began<-proc.time()[[3]]
  p<-withCallingHandlers(smashr::smash.gaus(y,sigma=sigma,v.est=FALSE,joint=FALSE,
    homoskedastic=FALSE,post.var=FALSE,reflect=TRUE,family='DaubExPhase',filter.number=4,
    ashparam=list(),jash=FALSE),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  stopifnot(!length(warnings),length(.GlobalEnv$.niv_smash_records)>0L)
  list(prediction=as.numeric(p),records=.GlobalEnv$.niv_smash_records,gridmult=64,
       options=expected,warnings=warnings,seconds=proc.time()[[3]]-began,valid=TRUE)
}

niv_check_smash <- function(z,y,sigma,D) {
  stopifnot(z$valid,!length(z$warnings),length(z$records)==D$J,
            identical(z$options,niv_options()),length(z$prediction)==D$n,
            all(is.finite(z$prediction)))
  theta<-numeric(length(D$aa));d<-drop(D$A%*%y)
  for(j in 0:(D$J-1L)) {
    r<-z$records[[j+1L]];ix<-j*D$N+seq_len(D$N);ix<-ix[D$native[ix]>0]
    niv_check_native(r$native);stopifnot(identical(r$options,z$options))
    noc_assert(r$wc,d[ix]);noc_assert(r$se^2,sigma^2*D$native[ix])
    g<-r$fit$fitted_g
    stopifnot(inherits(g,'normalmix'),all(g$mean==0),all(is.finite(c(g$pi,g$sd))),
              all(g$sd>=0),all(g$pi>=0),abs(sum(g$pi)-1)<=1e-8)
    pm<-niv_pm(r$wc,r$se^2,g$pi,g$sd);noc_assert(pm,ashr::get_pm(r$fit));theta[ix]<-pm
  }
  noc_assert(z$prediction,drop(D$synthesis%*%theta+D$scaling%*%y))
  invisible(TRUE)
}

nsm_geometry<-function(n,filter.number,family) {
  J<-as.integer(log2(n))
  stopifnot(n==2^J,J>=2L,filter.number==4,family=='DaubExPhase')
  e<-numeric(n);e[1]<-1
  w<-wavethresh::wd(e,filter.number=filter.number,family=family,type='station')
  # wd$D has exactly the row order of the original ndwt.mat()/wd.D operator.
  norms<-colSums(matrix(w$D^2,nrow=n,ncol=J))
  stopifnot(length(norms)==J,all(is.finite(norms)),all(norms>0))
  list(kind='native-periodic-row-norms-v02',n=n,J=J,norms=norms)
}

nsm_variance<-function(data.var,Wl,J,n) {
  stopifnot(identical(Wl$kind,'native-periodic-row-norms-v02'),Wl$n==n,Wl$J==J,
    length(data.var)==n,all(is.finite(data.var)),all(data.var>0),
    all(data.var==data.var[1]))
  # Native SMASH variances, deliberately not the reflected NOC covariance.
  rep(data.var[1]*Wl$norms,each=n)
}

nsm_rewrite<-function(fun,replacements,env) {
  counts<-integer(length(replacements))
  walk<-function(e) {
    for(i in seq_along(replacements))if(identical(e,replacements[[i]][[1]])) {
      counts[i]<<-counts[i]+1L;return(replacements[[i]][[2]])
    }
    if(is.call(e)) {
      parts<-as.list(e)
      for(k in seq_along(parts))if(!identical(unname(parts[k]),unname(alist(x=))))
        parts[k]<-list(walk(parts[[k]]))
      e<-as.call(parts)
    }
    e
  }
  body(fun)<-walk(body(fun));stopifnot(all(counts==1L));environment(fun)<-env;fun
}

nsm_functions<-function() {
  ns<-asNamespace('smashr');private<-new.env(parent=ns)
  private$nsm_geometry<-nsm_geometry;private$nsm_variance<-nsm_variance
  private$ndwt.mat<-function(...)stop('Dense wavelet geometry is forbidden in this adapter')
  private$mu.smooth<-nsm_rewrite(get('mu.smooth',ns),list(list(
    quote(apply((rep(1,n*J)%o%data.var)*Wl$W2,1,sum)),
    quote(nsm_variance(data.var,Wl,J,n)))),private)
  private$smash.gaus<-nsm_rewrite(get('smash.gaus',ns),list(list(
    quote(ndwt.mat(n,filter.number=filter.number,family=family)),
    quote(nsm_geometry(n,filter.number,family)))),private)
  adapter<-new.env(parent=.GlobalEnv)
  adapter$nsm_gaus<-function(x,sigma,v.est,joint,homoskedastic,post.var,reflect,
    family,filter.number,ashparam,jash) {
    stopifnot(length(sigma)==1L,is.finite(sigma),sigma>0,
      identical(v.est,FALSE),identical(joint,FALSE),identical(homoskedastic,FALSE),
      identical(post.var,FALSE),identical(reflect,TRUE),family=='DaubExPhase',
      filter.number==4,identical(ashparam,list()),identical(jash,FALSE))
    private$smash.gaus(x,sigma=sigma,v.est=v.est,joint=joint,
      homoskedastic=homoskedastic,post.var=post.var,reflect=reflect,
      family=family,filter.number=filter.number,ashparam=ashparam,jash=jash)
  }
  adapter$niv_smash<-nsm_rewrite(niv_smash,list(list(quote(smashr::smash.gaus),
    quote(nsm_gaus))),adapter)
  adapter$private<-private
  adapter
}

nsm_smash<-function(y,sigma)nsm_functions()$niv_smash(y,sigma)

pc_equal <- function(x, y, tol = 1e-12) {
  isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
}

pc_kernel_curve <- function(dat, h, ngrid) {
  KernSmooth::locpoly(dat$x, dat$y, drv = 0L, degree = 1L, kernel = "normal",
    bandwidth = h, gridsize = ngrid, range.x = c(0, 1), binned = FALSE, truncate = TRUE)
}

pc_ss_unpack <- function(v, nk) {
  b<-matrix(v,nk,4); M<-diag(b[,1])
  for(k in 1:3) {
    i<-seq_len(nk-k)
    M[cbind(i,i+k)]<-M[cbind(i+k,i)]<-b[i,k+1]
  }
  M
}

pc_ss_system <- function(reference, dat, cfg) {
  n<-length(dat$x); nk<-reference$fit$nk; tol<-cfg$ss_stable$matrix_tol
  stopifnot(nk==n+2L,pc_equal(reference$x,dat$x),pc_equal(reference$yin,dat$y),
    all(reference$w==1),reference$spar==cfg$ss_stable$probe_spar,
    is.finite(reference$ratio),reference$ratio>0,length(reference$auxM$Sigma)==4L*nk)
  er<-reference$iparms["errorI"]
  stopifnot(length(er)==1L,is.na(er)||er==0L)
  x<-(dat$x-reference$fit$min)/reference$fit$range
  X<-splines::splineDesign(reference$fit$knot,x,ord=4L)
  S<-pc_ss_unpack(reference$auxM$Sigma,nk)
  stopifnot(max(abs(crossprod(X)-pc_ss_unpack(reference$auxM$XWX,nk)))<=tol,
    max(abs(crossprod(X,dat$y)-reference$auxM$XWy))<=tol*(1+max(abs(dat$y))))
  # H maps node values to the minimum-native-penalty B-spline coefficients.
  q<-qr(t(X)); stopifnot(q$rank==n)
  Q<-qr.Q(q,complete=TRUE); T<-Q[,seq_len(n)]%*%solve(t(qr.R(q)))
  N<-Q[,(n+1L):nk,drop=FALSE]
  H<-T-N%*%solve(crossprod(N,S%*%N),crossprod(N,S%*%T))
  K<-crossprod(H,S%*%H)
  affine<-cbind(1,x)
  stopifnot(max(abs(X%*%H-diag(n)))<=tol,
    max(abs(K%*%affine))/(1+max(abs(K)))<=tol)
  # Preserve the affine null space explicitly instead of estimating its two zeros.
  V<-qr.Q(qr(affine),complete=TRUE)[,3:n,drop=FALSE]
  E<-eigen(crossprod(V,K%*%V),symmetric=TRUE)
  stopifnot(all(is.finite(E$values)),all(E$values>0))
  U<-V%*%E$vectors; d<-E$values; z<-drop(crossprod(U,dat$y))
  list(H=H,K=K,U=U,d=d,z=z,n=n,ratio=reference$ratio,
    affine=affine,reference=reference)
}

pc_ss_score <- function(system, spar) {
  lambda<-system$ratio*256^(3*spar-1)
  v<-system$d/(1+lambda*system$d)
  # Cancel lambda BEFORE division; no y - fitted or 1 - df/n cancellation.
  system$n*sum((v*system$z)^2)/sum(v)^2
}

pc_ss_search <- function(system,cfg) {
  s<-seq(cfg$ss$low,cfg$ss$high,length.out=cfg$ss_stable$search_points)
  g<-vapply(s,function(t) pc_ss_score(system,t),0.0)
  stopifnot(all(is.finite(g)),all(g>=0))
  j<-2:(length(s)-1L)
  j<-j[g[j]<=g[j-1L] & g[j]<=g[j+1L] & (g[j]<g[j-1L] | g[j]<g[j+1L])]
  intervals<-rbind(c(cfg$ss$low,cfg$ss$high),cbind(s[j-1L],s[j+1L]))
  refined<-lapply(seq_len(nrow(intervals)),function(i)
    stats::optimize(function(t) pc_ss_score(system,t),intervals[i,],tol=cfg$ss$tol))
  candidates<-data.frame(spar=c(s,vapply(refined,`[[`,0.0,"minimum")),
    score=c(g,vapply(refined,`[[`,0.0,"objective")))
  # Exact numerical ties select stronger smoothing. No truth-dependent tie tolerance.
  candidates<-candidates[order(-candidates$spar,method="radix"),];rownames(candidates)<-NULL
  selected<-which.min(candidates$score)
  list(grid=data.frame(spar=s,score=g),intervals=intervals,candidates=candidates,
    selected=selected,spar=candidates$spar[selected],score=candidates$score[selected])
}

pc_ss_solution <- function(system,spar,grid,dat) {
  lambda<-system$ratio*256^(3*spar-1)
  a<-lambda*system$d/(1+lambda*system$d)
  residual<-drop(system$U%*%(a*system$z))
  fitted<-dat$y-residual
  coefficients<-drop(system$H%*%fitted)
  f<-system$reference$fit;f$coef<-coefficients
  list(spar=spar,lambda=lambda,df=2+sum(1/(1+lambda*system$d)),
    residual_df=sum(a),gcv=pc_ss_score(system,spar),rss=sum((a*system$z)^2),
    fitted=fitted,residual=residual,coefficients=coefficients,fit=f,
    prediction=drop(stats::predict(f,grid)$y),
    spectral_stationarity=max(abs(a*system$z-lambda*system$d*system$z/(1+lambda*system$d))))
}

pc_ss_stable_fit <- function(dat,grid,cfg) {
  reference<-stats::smooth.spline(dat$x,dat$y,w=rep(1,length(dat$x)),all.knots=TRUE,
    spar=cfg$ss_stable$probe_spar,cv=FALSE,keep.data=TRUE,keep.stuff=TRUE,df.offset=0,penalty=1)
  system<-pc_ss_system(reference,dat,cfg)
  search<-pc_ss_search(system,cfg)
  solution<-pc_ss_solution(system,search$spar,grid,dat)
  out<-list(backend="native-penalty-stable-GCV-v02",reference=reference,search=search,
    solution=solution,prediction=solution$prediction,
    selected=list(criterion="stable_GCV",spar=solution$spar,lambda=solution$lambda,
      df=solution$df,residual_df=solution$residual_df,criterion_value=solution$gcv,
      boundary=min(abs(solution$spar-c(cfg$ss$low,cfg$ss$high)))<=cfg$ss_boundary_tol))
  pc_ss_stable_validate(out,dat,grid,cfg)
  out
}

pc_ss_stable_validate <- function(out,dat,grid,cfg) {
  stopifnot(identical(out$backend,"native-penalty-stable-GCV-v02"))
  system<-pc_ss_system(out$reference,dat,cfg)
  search<-pc_ss_search(system,cfg)
  solution<-pc_ss_solution(system,search$spar,grid,dat)
  stopifnot(pc_equal(out$search,search),pc_equal(out$solution,solution),
    pc_equal(out$prediction,solution$prediction),all(is.finite(solution$prediction)),
    all(is.finite(solution$coefficients)),solution$lambda>0,
    solution$df>=2-1e-6,solution$df<=length(dat$x)+1e-6,
    solution$residual_df>0,solution$gcv>=0,is.finite(solution$gcv),
    solution$spectral_stationarity<=cfg$ss_prediction_tol*(1+max(abs(dat$y))),
    max(abs(stats::predict(solution$fit,dat$x)$y-solution$fitted))<=
      cfg$ss_prediction_tol*(1+max(abs(solution$fitted))))
  expected<-list(criterion="stable_GCV",spar=solution$spar,lambda=solution$lambda,
    df=solution$df,residual_df=solution$residual_df,criterion_value=solution$gcv,
    boundary=min(abs(solution$spar-c(cfg$ss$low,cfg$ss$high)))<=cfg$ss_boundary_tol)
  stopifnot(pc_equal(out$selected,expected))
  # The probe checks equivalence where native arithmetic is well conditioned.
  probe<-pc_ss_solution(system,cfg$ss_stable$probe_spar,dat$x,dat)
  stopifnot(max(abs(probe$prediction-out$reference$y))<=
    cfg$ss_prediction_tol*(1+max(abs(out$reference$y))))
  invisible(TRUE)
}

nr_kernel <- function(d,h) {
  stopifnot(length(h)==1L,is.finite(h),h>0)
  grids<-c(4097L,8193L)
  fits<-lapply(grids,function(g)pc_kernel_curve(d,h,g))
  nr_raw(paste0('kernel_h_',format(h,digits=17)),fits)
  predictions<-lapply(seq_along(grids),function(j) {
    index<-1+d$x*(grids[j]-1L)
    stopifnot(all(index==as.integer(index)),pc_equal(fits[[j]]$x[index],d$x),
      all(is.finite(fits[[j]]$y)))
    fits[[j]]$y[index]
  })
  error<-max(abs(predictions[[1]]-predictions[[2]]))/(1+max(abs(predictions[[2]])))
  stopifnot(is.finite(error),error<=1e-4)
  list(prediction=predictions[[2]],bandwidth=h,resolution_error=error,grids=grids)
}

nr_dpill <- function(d) {
  h<-KernSmooth::dpill(d$x,d$y,blockmax=5,divisor=20,trim=.01,
    proptrun=.05,gridsize=401L,range.x=c(0,1),truncate=TRUE)
  nr_raw('kernel_dpill',h);stopifnot(length(h)==1L,is.finite(h),h>0);h
}

iv_gam_system <- function(f,d) {
  X<-predict(f,newdata=data.frame(x=d$x),type="lpmatrix")
  S<-matrix(0,ncol(X),ncol(X));sm<-f$smooth[[1]];ix<-sm$first.para:sm$last.para
  for(j in seq_along(sm$S))S[ix,ix]<-S[ix,ix]+f$sp[j]*sm$S[[j]]
  list(X=X,S=S,G=crossprod(X),rhs=drop(crossprod(X,d$y)),y=d$y)
}

ncv5_gradient <- function(f) {
  stopifnot(identical(f$method, 'REML'), identical(f$optimizer, c('outer', 'newton')),
    is.finite(f$sig2), f$sig2 > 0, identical(f$sig2, f$scale),
    !isTRUE(f$control$edge.correct))
  score <- as.numeric(f$gcv.ubre); epsilon <- f$control$newton$conv.tol
  stopifnot(is.finite(score), is.finite(epsilon), epsilon > 0)
  score.scale <- abs(log(f$sig2)) + abs(score)
  legacy.limit <- 5 * (abs(f$sig2) + abs(score)) * epsilon
  limit <- 5 * score.scale * epsilon; grad <- max(abs(f$outer.info$grad))
  hist <- f$outer.info$score.hist
  change <- if (length(hist) >= 2L) abs(tail(hist, 1L) - tail(hist, 2L)[1L]) else NA_real_
  list(rule='mgcv 1.9-4 REML/Newton score scale; versioned v05 correction',
    gradient=grad, relative_tolerance=epsilon, score_scale=score.scale,
    gradient_limit=limit, gradient_pass=grad <= limit,
    legacy_gradient_limit=legacy.limit, legacy_gradient_pass=grad <= legacy.limit,
    final_score_change=change, score_change_limit=score.scale * epsilon)
}

cp_gam_check <- function(f,d,k,method) {
  # Other optimizer/method configurations have no approved gate correction.
  if (!identical(f$method, 'REML') || !identical(f$optimizer, c('outer','newton')))
    return(ncv5_legacy_gam_check(f,d,k,method))
  gate <- ncv5_gradient(f)
  X<-predict(f,newdata=data.frame(x=d$x),type='lpmatrix');beta<-coef(f)
  S<-matrix(0,length(beta),length(beta));sm<-f$smooth[[1]]
  ix<-sm$first.para:sm$last.para
  for(j in seq_along(sm$S))S[ix,ix]<-S[ix,ix]+f$sp[j]*sm$S[[j]]
  g<-drop(crossprod(X,d$y-drop(X%*%beta))-S%*%beta)
  normal<-max(abs(g))/(1+max(abs(crossprod(X,d$y)),abs(S%*%beta)))
  ev<-eigen(f$outer.info$hess,symmetric=TRUE,only.values=TRUE)$values
  grad<-gate$gradient;tol<-gate$gradient_limit
  ck<-list(converged=isTRUE(f$converged),outer=f$outer.info$conv,
    normal_equation_residual=normal,edf=sum(f$edf),coefficient_dimension=length(beta),
    gradient=grad,native_gradient_tolerance=tol,hessian_eigenvalues=ev,
    relative_min_hessian=min(ev)/max(1,max(abs(ev))),scale=f$sig2,
    unsupported_columns=sum(colSums(abs(X))==0),sp=f$sp,
    basis=k,method=method,control=f$control,REML_gate_v05=gate)
  stopifnot(ck$converged,identical(ck$outer,'full convergence'),normal<=1e-8,
    all(is.finite(c(beta,f$sp,ev,grad,tol,f$sig2,ck$edf))),all(f$sp>0),f$sig2>0,
    grad<=tol,ck$relative_min_hessian>= -1e-8,ck$edf<=length(beta)+1e-8)
  ck
}

ncv6_previous_check <- cp_gam_check

ap4_payload <- function(sys,c,z) {
  hex<-function(x)I(sprintf('%a',as.vector(x)))
  list(n=nrow(sys$X),p=ncol(sys$X),X=hex(sys$X),S=hex(sys$S),rhs=hex(sys$rhs),
    y=hex(sys$y),hi=hex(z$coefficients),lo=hex(z$coefficient_correction),
    prediction=hex(z$prediction),multiplier=sprintf('%a',c),tolerance=sprintf('%a',1e-8))
}

ncv6_precision <- function(sys,c,method,phase,rejection) {
  stopifnot(method %in% c('PS','AP'),phase %in% c('native','fixed'))
  event<-list(method=method,phase=phase,basis=ncol(sys$X),multiplier=c,
    operation='shared two-part coefficient solve and exact residual certificate',
    reserved_before_execution=TRUE)
  audit<-getOption('noc.psap.precision.audit')
  if(!is.null(audit))audit$events[[length(audit$events)+1L]]<-event
  home<-getOption('noc.r.reference.raw')
  if(!is.null(home))cat(jsonlite::toJSON(event,auto_unbox=TRUE),'\n',
    file=file.path(home,'precision_ledger_v06.jsonl'),append=TRUE)
  z<-ap4_solve_cpp(sys$X,sys$S,sys$G,sys$rhs,c)
  cert<-ap4_certificate(sys,c,z)
  out<-list(prediction=as.numeric(z$prediction),coefficients=as.numeric(z$coefficients),
    coefficient_correction=as.numeric(z$coefficient_correction),
    normal_equation_residual=cert$normal_equation_residual_upper,multiplier=c,valid=TRUE,
    REML_reoptimized=FALSE,numerical_backend='shared PS/AP two-part solve; exact dyadic certificate v06',
    precision_certificate=cert,original_rejection=rejection,
    coefficient_representation=z$representation,penalty_representation=z$penalty_representation,
    ridge_added=z$ridge_added,penalty_clipping=z$penalty_clipping)
  nr_raw(paste0('PSAP_',method,'_',phase,'_K',ncol(sys$X),'_c',format(c,digits=17)),out)
  out
}

ncv6_verify_solution <- function(sys,c,z) {
  stopifnot(z$valid,identical(as.numeric(z$multiplier),as.numeric(c)),
    !z$REML_reoptimized,z$ridge_added==0,!z$penalty_clipping)
  cert<-ap4_certificate(sys,c,z)
  stopifnot(identical(cert,z$precision_certificate),
    identical(z$normal_equation_residual,cert$normal_equation_residual_upper))
  invisible(TRUE)
}

ncv6_native_check <- function(f,d,k,method,solution=NULL,allow_solve=FALSE) {
  stopifnot(identical(as.numeric(f$y),as.numeric(d$y)),
    identical(as.numeric(f$model$x),as.numeric(d$x)))
  original<-tryCatch(ncv6_previous_check(f,d,k,method),error=identity)
  if(!inherits(original,'error')) {
    stopifnot(is.null(solution))
    return(original)
  }
  stopifnot(method %in% c('PS','AP'),
    identical(conditionMessage(original),'normal <= 1e-08 is not TRUE'))
  # Every non-coefficient gate is checked on the unchanged native mgcv object.
  gate<-ncv5_gradient(f);beta<-coef(f)
  ev<-eigen(f$outer.info$hess,symmetric=TRUE,only.values=TRUE)$values
  stopifnot(isTRUE(f$converged),identical(f$outer.info$conv,'full convergence'),
    all(is.finite(c(beta,f$sp,ev,gate$gradient,gate$gradient_limit,f$sig2,sum(f$edf)))),
    all(f$sp>0),f$sig2>0,gate$gradient_pass,
    min(ev)/max(1,max(abs(ev)))>= -1e-8,sum(f$edf)<=length(beta)+1e-8)
  sys<-iv_gam_system(f,d)
  raw.g<-drop(crossprod(sys$X,d$y-drop(sys$X%*%beta))-sys$S%*%beta)
  raw.normal<-max(abs(raw.g))/(1+max(abs(sys$rhs),abs(sys$S%*%beta)))
  if(is.null(solution)) {
    stopifnot(allow_solve)
    solution<-ncv6_precision(sys,1,method,'native',conditionMessage(original))
  }
  ncv6_verify_solution(sys,1,solution)
  list(converged=isTRUE(f$converged),outer=f$outer.info$conv,
    normal_equation_residual=solution$normal_equation_residual,
    edf=sum(f$edf),coefficient_dimension=length(beta),gradient=gate$gradient,
    native_gradient_tolerance=gate$gradient_limit,hessian_eigenvalues=ev,
    relative_min_hessian=min(ev)/max(1,max(abs(ev))),scale=f$sig2,
    unsupported_columns=sum(colSums(abs(sys$X))==0),sp=f$sp,basis=k,method=method,
    control=f$control,REML_gate_v05=gate,native_raw_normal_equation_residual=raw.normal,
    native_raw_rejection=conditionMessage(original),coefficient_solution_v06=solution)
}

cp_gam_check <- function(f,d,k,method)ncv6_native_check(f,d,k,method,allow_solve=TRUE)

ncv6_anchor <- function(f,checks,d) {
  raw<-as.numeric(predict(f,newdata=data.frame(x=d$x),type='response'))
  z<-checks$coefficient_solution_v06
  if(is.null(z))return(list(prediction=raw,fit=f,checks=checks))
  list(prediction=z$prediction,fit=f,checks=checks,coefficient_solution=z,
    native_raw_prediction=raw,REML_reoptimized=FALSE)
}

cp_gam <- function(d,k,method) {
  # Same mgcv call as frozen cp_gam; raw f is never edited by refinement.
  dat<-data.frame(x=d$x,y=d$y)
  form<-if(method=='PS')y~s(x,bs='ps',k=k,m=c(2,2)) else
    y~s(x,bs='ad',k=k,m=5,xt=list(bs='ps'))
  environment(form)<-environment();s<-mgcv::s
  fit<-mgcv::gam(form,data=dat,family=gaussian(),method='REML')
  ck<-cp_gam_check(fit,d,k,method)
  ncv6_anchor(fit,ck,d)
}

ncv6_validate_anchor <- function(a,d,K,method) {
  ck<-ncv6_native_check(a$fit,d,K,method,a$coefficient_solution,allow_solve=FALSE)
  if(is.null(a$coefficient_solution)) {
    stopifnot(identical(a$prediction,as.numeric(predict(a$fit,newdata=data.frame(x=d$x),type='response'))))
  } else {
    stopifnot(identical(ck,a$checks),identical(a$prediction,a$coefficient_solution$prediction),
      identical(a$native_raw_prediction,as.numeric(predict(a$fit,newdata=data.frame(x=d$x),type='response'))),
      !a$REML_reoptimized)
  }
  invisible(TRUE)
}

mc3_policy<-list(version='native-REML-saved-fit-v03',basis=96L,normal_limit=1e-8,
  hessian_flag_limit=-1e-8,hessian_is_acceptance_veto=FALSE,
  native_gradient_and_score_change=TRUE,new_REML_fits=0L,new_noise_draws=0L)

mc3_check<-function(f,d,k,method,solution=NULL,allow_solve=FALSE) {
  stopifnot(method%in%c('PS','AP'),k==96L,length(coef(f))==96L,
    length(f$sp)==if(method=='PS')1L else 5L,
    identical(as.numeric(f$y),as.numeric(d$y)),identical(as.numeric(f$model$x),as.numeric(d$x)),
    identical(f$family$family,'gaussian'),identical(f$family$link,'identity'))
  g<-ncv5_gradient(f);beta<-coef(f);H<-f$outer.info$hess
  stopifnot(is.matrix(H),nrow(H)==length(f$sp)+1L,ncol(H)==nrow(H),all(is.finite(H)))
  ev<-eigen(H,symmetric=TRUE,only.values=TRUE)$values;rel<-min(ev)/max(1,max(abs(ev)))
  stopifnot(isTRUE(f$converged),identical(f$outer.info$conv,'full convergence'),
    all(is.finite(c(beta,f$sp,ev,f$sig2,sum(f$edf),g$gradient,g$gradient_limit,
      g$final_score_change,g$score_change_limit))),all(f$sp>0),f$sig2>0,
    g$gradient_pass,g$final_score_change<=g$score_change_limit,sum(f$edf)<=length(beta)+1e-8)
  sys<-iv_gam_system(f,d)
  r<-drop(crossprod(sys$X,d$y-drop(sys$X%*%beta))-sys$S%*%beta)
  normal<-max(abs(r))/(1+max(abs(sys$rhs),abs(sys$S%*%beta)))
  if(is.null(solution)&&normal>1e-8) {
    if(!allow_solve)stop('Saved coefficients fail 1e-8; fixed-penalty correction required')
    solution<-ncv6_precision(sys,1,method,'native','v03 saved coefficient residual exceeds 1e-8')
  }
  if(!is.null(solution))ncv6_verify_solution(sys,1,solution)
  accepted_normal<-if(is.null(solution))normal else solution$normal_equation_residual
  stopifnot(is.finite(accepted_normal),accepted_normal<=1e-8)
  list(converged=TRUE,outer=f$outer.info$conv,normal_equation_residual=accepted_normal,
    native_raw_normal_equation_residual=normal,edf=sum(f$edf),coefficient_dimension=length(beta),
    gradient=g$gradient,native_gradient_tolerance=g$gradient_limit,REML_gate_v05=g,
    hessian_eigenvalues=ev,relative_min_hessian=rel,hessian_flag=rel< -1e-8,
    hessian_is_acceptance_veto=FALSE,scale=f$sig2,sp=f$sp,basis=k,method=method,
    control=f$control,coefficient_solution_v06=solution,acceptance_policy=mc3_policy)
}

mc3_validate_anchor<-function(a,d,K,method) {
  ck<-mc3_check(a$fit,d,K,method,a$coefficient_solution,FALSE)
  stopifnot(identical(ck,a$checks),identical(a$prediction,ncv6_anchor(a$fit,ck,d)$prediction),
    all(is.finite(a$prediction)),length(a$prediction)==length(d$y))
  invisible(TRUE)
}

cvp_record <- function(m,d,p,t,value=t$value,native=rep(FALSE,ncol(p)),checks=list()) {
  p<-as.matrix(p);stopifnot(nrow(p)==length(d$y),ncol(p)==nrow(t),
    length(value)==ncol(p),length(native)==ncol(p),all(is.finite(p)),
    !anyDuplicated(t$candidate))
  t$value<-value;t$native<-native
  list(method=m,id=d$id,prediction=p,table=t,n=length(d$y),checks=checks,
    sigma_hat=d$sigma_hat,valid=TRUE,truth_read=FALSE,heldout_read=FALSE)
}

cvp_ss <- function(d,t) {
  ref<-stats::smooth.spline(d$x,d$y,w=rep(1,length(d$x)),all.knots=TRUE,
    spar=.5,cv=FALSE,keep.data=TRUE,keep.stuff=TRUE,df.offset=0,penalty=1)
  sys<-pc_ss_system(ref,d,cp_cfg)
  probe<-pc_ss_solution(sys,.5,d$x,d)
  noc_assert(probe$prediction,ref$y)
  search<-pc_ss_search(sys,cp_cfg);v<-t$value
  v[is.na(v)]<-search$spar
  fits<-lapply(v,function(s) {
    z<-pc_ss_solution(sys,s,d$x,d)
    stopifnot(all(is.finite(c(z$prediction,z$df,z$gcv))),
      z$df>=2-1e-6,z$df<=length(d$x)+1e-6,z$residual_df>0,
      max(abs(z$prediction-z$fitted))<=1e-8*(1+max(abs(z$fitted))))
    z
  })
  cvp_record('SS',d,vapply(fits,function(a)a$prediction,numeric(length(d$y))),t,v,
    t$candidate=='native_GCV',list(df=vapply(fits,function(a)a$df,0.0),
      actual_lambda=vapply(fits,function(a)a$lambda,0.0)))
}

n3r5_grids<-as.integer(2^(12:20)+1)

n3r5_nodes<-function(d,h,M) {
  grid<-seq(0,1,length.out=M);delta<-1/(M-1L)
  b<-KernSmooth:::rlbin(d$x,d$y,grid,TRUE)
  occupied<-which(b$xcounts>0)
  location<-d$x/delta+1;lower<-floor(location);upper<-ceiling(location)
  stopifnot(all(lower>=1),all(upper<=M))
  wanted<-sort(unique(as.integer(c(lower,upper))))
  z<-n3r5_binned_nodes_cpp(as.integer(occupied),b$xcounts[occupied],b$ycounts[occupied],
    wanted,delta,h,as.integer(min(M-1L,floor(4*h/delta))))
  lo<-z$prediction[match(lower,wanted)];hi<-z$prediction[match(upper,wanted)]
  fraction<-location-lower
  prediction<-lo+fraction*(hi-lo)
  stopifnot(all(is.finite(prediction)),all(z$normal_residual<=1e-8))
  list(prediction=prediction,nodes=wanted,node_prediction=z$prediction,
    node_normal_residual=z$normal_residual,node_rank=z$rank,node_slope=z$slope,grid=M)
}

n3r5_pair<-function(d,h,candidate) {
  previous<-NULL;history<-list()
  for(M in n3r5_grids) {
    nr_raw('kernel_consumed_nodes',list(candidate=candidate,bandwidth=h,grid=M,trace=history))
    current<-n3r5_nodes(d,h,M)
    if(!is.null(previous)) {
      error<-max(abs(previous$prediction-current$prediction))/(1+max(abs(current$prediction)))
      history[[length(history)+1L]]<-list(coarse_grid=previous$grid,fine_grid=M,error=error)
      if(is.finite(error)&&error<=1e-4)return(list(prediction=current$prediction,
        coarse_prediction=previous$prediction,resolution_error=error,trace=history,
        consumed_coarse=previous,consumed_fine=current))
    }
    previous<-current
  }
  nr_raw('kernel_consumed_nodes',list(candidate=candidate,bandwidth=h,grid=M,trace=history))
  stop(sprintf('Consumed-node kernel resolution failed: %s h=%.17g final error=%.17g > 1e-4',candidate,h,error))
}

n3r5_kernel<-function(d,t) {
  h<-t$value;h[is.na(h)]<-nr_dpill(d)
  stopifnot(all(is.finite(h)),all(h>0))
  fits<-lapply(seq_along(h),function(j)n3r5_pair(d,h[j],t$candidate[j]))
  z<-cvp_record('KERNEL',d,vapply(fits,`[[`,numeric(length(d$y)),'prediction'),t,h,
    t$candidate=='native_dpill',list(resolution_error=vapply(fits,`[[`,0.0,'resolution_error'),
      coarse_prediction=vapply(fits,`[[`,numeric(length(d$y)),'coarse_prediction'),
      resolution_trace=lapply(fits,`[[`,'trace'),consumed_coarse=lapply(fits,`[[`,'consumed_coarse'),
      consumed_fine=lapply(fits,`[[`,'consumed_fine'),grids=n3r5_grids,
      evaluation_scope='all nodes consumed by original linear readout; unused nodes neither computed nor imputed',
      objective='unchanged KernSmooth linear bins, tau=4 Gaussian degree-one least squares',
      readout='same public-grid linear interpolation',bandwidth_candidates_changed=FALSE,tolerance_changed=FALSE))
  n3r5_kernel_check(z,d);z
}

n3r5_kernel_check<-function(g,d) {
  a<-g$checks;stopifnot(g$valid,g$n==length(d$y),!g$truth_read,!g$heldout_read,
    !a$bandwidth_candidates_changed,!a$tolerance_changed,identical(a$grids,n3r5_grids))
  for(j in seq_len(ncol(g$prediction))) {
    coarse<-a$consumed_coarse[[j]];fine<-a$consumed_fine[[j]]
    for(z in list(coarse,fine)) {
      loc<-d$x*(z$grid-1L)+1;wanted<-sort(unique(as.integer(c(floor(loc),ceiling(loc)))))
      stopifnot(identical(wanted,z$nodes),all(is.finite(z$node_prediction)),
        all(is.finite(z$node_normal_residual)),all(z$node_normal_residual<=1e-8),all(z$node_rank%in%1:2))
      p<-stats::approx((z$nodes-1)/(z$grid-1),z$node_prediction,xout=d$x,ties='ordered',rule=2)$y
      stopifnot(max(abs(p-z$prediction))<=1e-12*(1+max(abs(p))))
      bins<-KernSmooth:::rlbin(d$x,d$y,seq(0,1,length.out=z$grid),TRUE)
      occupied<-which(bins$xcounts>0);delta<-1/(z$grid-1L)
      cert<-n3r5_binned_certificate_cpp(as.integer(occupied),bins$xcounts[occupied],
        bins$ycounts[occupied],z$nodes,z$node_prediction,z$node_slope,delta,g$table$value[j],
        as.integer(min(z$grid-1L,floor(4*g$table$value[j]/delta))))
      stopifnot(all(is.finite(cert)),all(cert<=1e-8),identical(cert,z$node_normal_residual))
    }
    stopifnot(identical(fine$prediction,g$prediction[,j]),identical(coarse$prediction,a$coarse_prediction[,j]))
    error<-max(abs(coarse$prediction-fine$prediction))/(1+max(abs(fine$prediction)))
    trace<-a$resolution_trace[[j]]
    stopifnot(identical(error,a$resolution_error[j]),error<=1e-4,
      identical(error,tail(trace,1)[[1]]$error))
    for(k in seq_along(trace))stopifnot(trace[[k]]$coarse_grid==n3r5_grids[k],
      trace[[k]]$fine_grid==n3r5_grids[k+1L])
    if(length(trace)>1L)stopifnot(all(vapply(head(trace,-1),`[[`,0.0,'error')>1e-4))
  }
  invisible(TRUE)
}

