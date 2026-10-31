# Original CIS and three comparison estimators; refitting remains an internal step.
# No data-generation, repeated contamination or cross-validation is run here.
.cis_mve_nsamp <- 1000L # Internal MVE subsets; not contamination replications.
.cis_methods <- c('MLE','GammaRidge','RefitRidgeHKB','GammaCIS')
.cis_colors <- c('#444444','#B57916','#8261A8','#1764A0')

# ---------------------- Numerical and fitting functions ----------------
cis_make_design <- function(z, center=NULL, spread=NULL) {
  z <- as.matrix(z); storage.mode(z) <- 'double'
  if (is.null(center)) center <- colMeans(z)
  if (is.null(spread)) spread <- apply(z,2,sd)
  if (any(!is.finite(spread) | spread<=0)) stop('Constant/nonfinite predictor')
  X <- cbind('(Intercept)'=1,sweep(sweep(z,2,center,'-'),2,spread,'/'))
  list(X=X,center=center,spread=spread)
}
cis_positive_mean <- function(X,b) {
  mu <- as.vector(exp(X %*% b))
  if (any(!is.finite(mu) | mu<=0)) stop('Nonfinite/nonpositive prediction')
  mu
}
cis_gamma_fit <- function(X,y) {
  if (any(!is.finite(y) | y<=0)) stop('Gamma requires strictly positive y')
  if (qr(X)$rank<ncol(X) || nrow(X)<=ncol(X)) stop('Insufficient full-rank data')
  fit <- glm.fit(x=X,y=pmax(y,1e-8),family=Gamma(link='log'),
                 control=glm.control(maxit=200,epsilon=1e-9))
  class(fit) <- c('glm','lm')
  fit$x <- X
  b <- coef(fit)
  if (!isTRUE(fit$converged) || any(!is.finite(b))) stop('Gamma fit did not converge')
  mu <- cis_positive_mean(X,b)
  phi <- sum(((y-mu)/mu)^2)/(nrow(X)-ncol(X))
  if (!is.finite(phi) || phi<=0) stop('Invalid Pearson dispersion')
  rdev <- sign(y-mu)*sqrt(pmax(0,2*(y/mu-1-log(y/mu))))
  list(beta=b,mu=mu,phi=phi,rdev=rdev,fit=fit)
}
# Solve the original quadratic transformation by scaled augmented QR.
# This changes the numerical implementation, not the estimator.
# Partial-pivoted elimination uses MPFR numbers throughout.
cis_mp_solve <- function(A,rhs) {
 q<-length(rhs);U<-A;v<-rhs
 for(j in seq_len(q)) {
  ix<-j:q;vabs<-abs(U[ix,j]);pv<-ix[which(vabs==max(vabs))[1L]]
  if(U[pv,j]==0)stop("Zero MPFR pivot at this precision")
  if(pv!=j) {tmp<-U[j,];U[j,]<-U[pv,];U[pv,]<-tmp;tmp<-v[j];v[j]<-v[pv];v[pv]<-tmp}
  if(j<q) for(i in (j+1L):q) {
   a<-U[i,j]/U[j,j];U[i,j:q]<-U[i,j:q]-a*U[j,j:q];v[i]<-v[i]-a*v[j]
  }
 }
 ans<-v*0
 for(i in q:1L) {
  z<-if(i<q)sum(U[i,(i+1L):q]*ans[(i+1L):q])else 0
  ans[i]<-(v[i]-z)/U[i,i]
 }
 ans
}
cis_mp_shrink <- function(X,mu,b,k,precisions=c(256L,512L,1024L,2048L,4096L,8192L)) {
 if(!requireNamespace("Rmpfr",quietly=TRUE))stop("MPFR recovery requires Rmpfr")
 previous<-NULL;audit<-list();q<-ncol(X)
 for(bits in precisions) {
  one<-tryCatch({
   # Rebuild from unchanged X and mu, never from a rounded double D.
   xm<-Rmpfr::mpfr(X,bits);mm<-Rmpfr::mpfr(mu,bits)
   D<-Rmpfr::mpfr(matrix(0,q,q),bits)
   for(i in seq_len(q))for(j in i:q) {
    D[i,j]<-sum(xm[,i]*mm*xm[,j]);D[j,i]<-D[i,j]
   }
   bm<-Rmpfr::mpfr(b,bits);km<-Rmpfr::mpfr(k,bits);A<-D
   for(j in seq_len(q))A[j,j]<-A[j,j]+km
   rhs<-bm*0
   for(i in seq_len(q))rhs[i]<-sum(D[i,]*bm)
   z<-cis_mp_solve(A,rhs);res<-rhs*0
   for(i in seq_len(q))res[i]<-sum(A[i,]*z)-rhs[i]
   residual<-max(abs(res))/(max(abs(A))*q*max(abs(z))+max(abs(rhs)))
   change<-if(is.null(previous))Inf else max(abs(z-Rmpfr::mpfr(previous,bits)))/(1+max(abs(z)))
   list(z=z,residual=Rmpfr::asNumeric(residual),change=Rmpfr::asNumeric(change))
  },error=function(e)e)
  if(inherits(one,"error")) {
   audit[[length(audit)+1L]]<-data.frame(Bits=bits,Status=conditionMessage(one),RelativeChange=NA_real_,RelativeResidual=NA_real_)
   next # A zero pivot at one precision does not terminate escalation.
  }
  audit[[length(audit)+1L]]<-data.frame(Bits=bits,Status="Solved",RelativeChange=one$change,RelativeResidual=one$residual)
  if(!is.null(previous)&&is.finite(one$change)&&one$change<1e-10&&one$residual<1e-30) {
   z<-Rmpfr::asNumeric(one$z);if(any(!is.finite(z)))stop("Nonfinite converted MPFR result")
   attr(z,"numerics")<-list(Solver="MPFR",Bits=bits,RelativeChange=one$change,
    RelativeResidual=one$residual,Attempts=do.call(rbind,audit));return(z)
  }
  previous<-one$z
 }
 stop(paste("MPFR recovery did not pass two-precision agreement/residual checks through",max(precisions),"bits:",paste(vapply(audit,function(a)paste(a$Bits,a$Status,a$RelativeChange,a$RelativeResidual),character(1)),collapse="; ")))
}
cis_shrink <- function(X,mu,b,k,force_mpfr=FALSE) {
 if(any(!is.finite(X))||any(!is.finite(mu))||any(mu<=0)||any(!is.finite(b))||!is.finite(k)||k<=0)
  stop("Invalid shrinkage inputs")
 fast<-if(force_mpfr)NULL else tryCatch({
  B<-X*sqrt(mu);sc<-max(abs(B),sqrt(k));A<-rbind(B/sc,sqrt(k)/sc*diag(ncol(X)))
  rhs<-c(as.vector((B/sc)%*%b),rep(0,ncol(X)))
  sv<-svd(A,nu=0,nv=0)$d
  if(min(sv)/max(sv)<1e-10)stop("Augmented system requires higher precision")
  qa<-qr(A,LAPACK=FALSE,tol=1e-12)
  if(qa$rank<ncol(X))stop("Numerically deficient augmented QR")
  z<-as.vector(qr.coef(qa,rhs));if(any(!is.finite(z)))stop("Nonfinite QR result")
  attr(z,"numerics")<-list(Solver="QR",Bits=53L,RelativeChange=NA_real_,RelativeResidual=NA_real_,Attempts=NULL)
  z
 },error=function(e)NULL)
 if(!is.null(fast))return(fast)
 cis_mp_shrink(X,mu,b,k)
}
cis_screen_refit <- function(X,y,mle) {
  Z <- X[,-1,drop=FALSE]
  S <- MASS::cov.rob(Z,method='mve',nsamp=.cis_mve_nsamp)$cov
  rd2 <- mahalanobis(Z,apply(Z,2,median),S)
  wx <- as.integer(rd2<=qchisq(.975,ncol(Z)))
  wr <- as.integer(abs(mle$rdev)<=2.5)
  raw_w <- wx*wr; fallback <- sum(raw_w)<=ncol(X)+5
  w <- if (fallback) rep(1L,length(y)) else raw_w
  take <- which(w==1L)
  ref <- if (fallback) mle else cis_gamma_fit(X[take,,drop=FALSE],y[take])
  # Evaluate means only on retained rows, avoiding 0 * Inf.
  mu <- cis_positive_mean(X[take,,drop=FALSE],ref$beta)
  list(ref=ref,take=take,mu=mu,w=w,raw_w=raw_w,wx=wx,wr=wr,
       rd2=rd2,fallback=fallback)
}
cis_fit_methods <- function(X,y) {
  q <- ncol(X); out <- setNames(vector('list',length(.cis_methods)),.cis_methods)
  errors <- setNames(rep('',length(.cis_methods)),.cis_methods); meta <- list()
  attempt <- function(label,fun) {
    tryCatch(fun(),error=function(e){errors[label]<<-conditionMessage(e);NULL})
  }
  ml <- attempt('MLE',function() cis_gamma_fit(X,y))
  if (is.null(ml)) {errors[]<-errors['MLE'];return(list(beta=out,errors=errors,meta=meta))}
  out$MLE <- ml$beta
  kh <- q*ml$phi/sum(ml$beta^2)
  out['GammaRidge'] <- list(attempt('GammaRidge',function() cis_shrink(X,ml$mu,ml$beta,kh)))
  sc_error <- ''
  sc <- tryCatch(cis_screen_refit(X,y,ml),error=function(e){sc_error <<- conditionMessage(e);NULL})
  meta <- list(mle=ml,k_HKB=kh)
  if (is.null(sc)) {
    errors[c('RefitRidgeHKB','GammaCIS')]<-sc_error
  } else {
    Xr <- X[sc$take,,drop=FALSE]; b <- sc$ref$beta
    C <- median(colSums((Xr*sqrt(sc$mu))^2))
    kf <- q*sc$ref$phi/sum(b^2)
    kc <- C*sc$ref$phi/(sum(b^2)+1e-8)
    out['RefitRidgeHKB'] <- list(attempt('RefitRidgeHKB',function() cis_shrink(Xr,sc$mu,b,kf)))
    out['GammaCIS'] <- list(attempt('GammaCIS',function() cis_shrink(Xr,sc$mu,b,kc)))
    meta <- c(meta,list(screen=sc,CIS=C,k_CIS=kc,k_refit_HKB=kf))
  }
  meta$numerics<-lapply(out,function(b)if(is.null(b))NULL else attr(b,"numerics"))
  list(beta=out,errors=errors,meta=meta)
}
# Transform standardized coefficients back to original predictor units.
cis_original_beta <- function(b,center,spread) {
  c(b[1]-sum(b[-1]*center/spread),b[-1]/spread)
}
cis_gamma_unit_deviance <- function(y,mu) 2*(y/mu-1-log(y/mu))

# Same unpenalized Gamma MLE; QR-coordinate likelihood fallback if IRLS fails.
# This is a numerical recovery, not a robust or penalized replacement.
cis_gamma_fit_irls<-cis_gamma_fit
cis_gamma_fit<-function(X,y){
 initial<-tryCatch(cis_gamma_fit_irls(X,y),error=function(e)e)
 if(!inherits(initial,'error'))return(initial)
 if(any(y<=0)||qr(X)$rank<ncol(X))stop(conditionMessage(initial))
 qrX<-qr(X);Q<-qr.Q(qrX);R<-qr.R(qrX);n<-nrow(X)
 b0<-c(log(mean(y)),rep(0,ncol(X)-1));g0<-as.vector(R%*%b0)
 fn<-function(g){eta<-as.vector(Q%*%g);t<-log(y)-eta;if(any(t>700))return(Inf);mean(eta+exp(t))}
 gr<-function(g){eta<-as.vector(Q%*%g);as.vector(crossprod(Q,1-exp(log(y)-eta)))/n}
 op<-optim(g0,fn,gr,method='BFGS',control=list(maxit=3000,reltol=1e-14))
 g<-op$par
 # Newton refinement in orthogonal coordinates with likelihood backtracking.
 for(it in 1:50){
  eta<-as.vector(Q%*%g);ratio<-exp(log(y)-eta);score<-as.vector(crossprod(Q,ratio-1))
  if(max(abs(crossprod(X,ratio-1)))<1e-4)break
  H<-crossprod(Q*sqrt(ratio));step<-tryCatch(solve(H,score),error=function(e)NULL)
  if(is.null(step))break
  rate<-1;old<-fn(g)
  while(rate>2^-30 && (!is.finite(fn(g+rate*step)) || fn(g+rate*step)>old+1e-13))rate<-rate/2
  g<-g+rate*step
 }
 b<-as.vector(backsolve(R,g));names(b)<-colnames(X);mu<-cis_positive_mean(X,b)
 if(max(abs(crossprod(X,y/mu-1)))>1e-4)stop('MLE recovery score check failed')
 phi<-sum(((y-mu)/mu)^2)/(n-ncol(X))
 if(!is.finite(phi)||phi<=0)stop('Invalid recovered dispersion')
 list(beta=b,mu=mu,phi=phi,rdev=sign(y-mu)*sqrt(pmax(0,2*(y/mu-1-log(y/mu)))),fit=NULL)
}

# Consistent formatted Excel tables and plain-text reproducibility metadata.
cis_export_tables <- function(tabs,path) {
 wb <- openxlsx::createWorkbook()
 for(nm in names(tabs)) {
  openxlsx::addWorksheet(wb,nm)
  if(nrow(tabs[[nm]])) openxlsx::writeDataTable(wb,nm,tabs[[nm]],tableStyle="TableStyleMedium2")
  else openxlsx::writeData(wb,nm,tabs[[nm]])
  openxlsx::freezePane(wb,nm,firstRow=TRUE)
  openxlsx::setColWidths(wb,nm,seq_len(ncol(tabs[[nm]])),"auto")
 }
 openxlsx::saveWorkbook(wb,path,overwrite=TRUE)
}


# Public matrix API. X must include its unit intercept and have fixed scaling.
cis_fit <- function(X,y,screen_seed=NULL) {
 X<-as.matrix(X);storage.mode(X)<-"double";y<-as.numeric(y)
 if(nrow(X)!=length(y)||ncol(X)<2||any(!is.finite(X))||any(!is.finite(y))||any(y<=0))stop("Invalid design/response")
 if(any(X[,1]!=1))stop("The first design column must be the unit intercept")
 if(!is.null(screen_seed)) {
  had<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had)old<-get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(had)assign(".Random.seed",old,envir=.GlobalEnv)else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv))
  set.seed(screen_seed)
 }
 cis_fit_methods(X,y)
}
