cis_simulate <- function(run_mode="full",B_override=NA_integer_,output_dir=file.path(getwd(),"CIS_Simulation"),setting_ids=NULL) {
for(pkg in c("MASS","openxlsx","Rmpfr"))if(!requireNamespace(pkg,quietly=TRUE))stop(paste("Install",pkg))
# -------------------------- USER SETTINGS --------------------------
RUN_MODE <- run_mode
B_OVERRIDE <- B_override
OUTPUT_DIR <- output_dir
if(!RUN_MODE %in% c("pilot","full")) stop("Invalid RUN_MODE")
if(RUN_MODE=="pilot") {
 n_set<-50L;p_set<-c(4L,8L);rho_set<-.95;phi_set<-.5;B<-50L
} else {
 n_set<-c(50L,100L);p_set<-c(4L,8L);rho_set<-c(0,.90,.95,.99)
 phi_set<-c(.25,.50,1);B<-500L
}
if(!is.na(B_OVERRIDE)) B<-as.integer(B_OVERRIDE)
stopifnot(B>=2)
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
profiles<-c("EqualPositive","MixedSigns")
# Explicit ordering retained in the saved settings; replication seeds depend on it.
base<-expand.grid(n=n_set,p=p_set,rho=rho_set,phi=phi_set,
 BetaProfile=profiles,stringsAsFactors=FALSE)
conditions<-data.frame(Scenario=c("Clean","YOutliers","BadLeverage","YOutliers","BadLeverage"),
 Rate=c(0,.1,.1,.2,.2))
settings<-do.call(rbind,lapply(seq_len(nrow(base)),function(i)cbind(base[rep(i,nrow(conditions)),,drop=FALSE],conditions)))
settings$SettingID<-seq_len(nrow(settings));rownames(settings)<-NULL
saveRDS(settings,file.path(OUTPUT_DIR,"Settings.rds"))
comparison<-performance<-paired<-screening<-errors<-list();numerical<-list()
if(!is.null(setting_ids)) {
 if(any(!setting_ids%in%settings$SettingID))stop("Unknown setting_ids")
 settings<-settings[settings$SettingID%in%setting_ids,,drop=FALSE]
}
for(s in seq_len(nrow(settings))) {
 st<-settings[s,];n<-st$n;p<-st$p;q<-p+1L
 beta<-c(1,if(st$BetaProfile=="EqualPositive")rep(.5,p) else rep(c(.6,-.4,.3,-.2),length.out=p))
 estimates<-array(NA_real_,c(B,length(.cis_methods),q),dimnames=list(NULL,.cis_methods,c("Intercept",paste0("x",1:p))))
 loss<-matrix(NA_real_,B,length(.cis_methods),dimnames=list(NULL,.cis_methods))
 scr<-data.frame(Replication=seq_len(B),UsedN=NA_real_,CIS=NA_real_,kCIS=NA_real_,kRHKB=NA_real_,Fallback=NA,
 Detection=NA_real_,FalseExclusion=NA_real_)
 err<-list()
 for(b in seq_len(B)) {
  seed<-123L+(st$SettingID-1L)*100000L+b;set.seed(seed)
  dat<-tryCatch({
   Sigma<-matrix(st$rho,p,p);diag(Sigma)<-1
   X<-cbind(1,scale(MASS::mvrnorm(n,rep(0,p),Sigma)))
   mu<-as.vector(exp(X%*%beta));y<-rgamma(n,shape=1/st$phi,scale=mu*st$phi)
   ids<-if(st$Rate>0)sample.int(n,ceiling(n*st$Rate))else integer()
   y[ids]<-y[ids]*8
   if(st$Scenario=="BadLeverage") X[ids,-1]<-X[ids,-1,drop=FALSE]*50
   list(X=X,y=y,ids=ids)
  },error=function(e)e)
  if(inherits(dat,"error")) {
   err[[length(err)+1L]]<-data.frame(Replication=b,Method="Generation",Message=conditionMessage(dat));next
  }
  obj<-tryCatch(suppressWarnings(cis_fit_methods(dat$X,dat$y)),error=function(e)e)
  if(inherits(obj,"error")) {
   err[[length(err)+1L]]<-data.frame(Replication=b,Method="All",Message=conditionMessage(obj));next
  }
  for(m in .cis_methods) {
   if(!is.null(obj$beta[[m]])) {
    estimates[b,m,]<-obj$beta[[m]];loss[b,m]<-mean((obj$beta[[m]]-beta)^2)
   } else err[[length(err)+1L]]<-data.frame(Replication=b,Method=m,Message=obj$errors[m])
  }
  for(m in .cis_methods) {
   nu<-obj$meta$numerics[[m]]
   if(!is.null(nu)) numerical[[length(numerical)+1L]]<-data.frame(SettingID=st$SettingID,Replication=b,Method=m,
    Solver=nu$Solver,Bits=nu$Bits,RelativeChange=nu$RelativeChange,RelativeResidual=nu$RelativeResidual)
  }
  sc<-obj$meta$screen
  if(!is.null(sc)) {
   scr$UsedN[b]<-sum(sc$w);scr$CIS[b]<-obj$meta$CIS;scr$kCIS[b]<-obj$meta$k_CIS
   scr$kRHKB[b]<-obj$meta$k_refit_HKB;scr$Fallback[b]<-sc$fallback
   if(length(dat$ids))scr$Detection[b]<-mean(sc$raw_w[dat$ids]==0)
   clean<-setdiff(seq_len(n),dat$ids);scr$FalseExclusion[b]<-mean(sc$raw_w[clean]==0)
  }
 }
 common<-apply(is.finite(loss),1,all)
 avg<-if(any(common))colMeans(loss[common,,drop=FALSE])else setNames(rep(NA_real_,length(.cis_methods)),.cis_methods)
 comparison[[s]]<-cbind(st,CommonSuccessful=sum(common),CommonFailRate=mean(!common),
  as.data.frame(as.list(setNames(avg,paste0("MSE_",.cis_methods)))),
  Best=if(all(is.na(avg)))NA_character_ else names(which.min(avg)))
 pp<-list();dd<-list()
 for(m in .cis_methods) {
  ok<-is.finite(loss[,m]);z<-matrix(estimates[ok,m,],ncol=q)
  bias2<-variance<-NA_real_
  if(sum(ok)) {bm<-colMeans(z);bias2<-mean((bm-beta)^2);variance<-mean(colMeans(sweep(z,2,bm,"-")^2))}
  pp[[m]]<-cbind(st,Estimator=m,Successful=sum(ok),FailRate=mean(!ok),Bias2=bias2,Variance=variance,
   MSE=if(any(ok))mean(loss[ok,m])else NA_real_)
  if(m!="GammaCIS") {
   ix<-is.finite(loss[,m])&is.finite(loss[,"GammaCIS"]);d<-loss[ix,m]-loss[ix,"GammaCIS"]
   se<-if(length(d)>1)sd(d)/sqrt(length(d))else NA_real_;dm<-if(length(d))mean(d)else NA_real_
   dd[[m]]<-cbind(st,Comparator=m,PairSuccessful=length(d),MeanDifference=dm,MCSE=se,Lower=dm-1.96*se,Upper=dm+1.96*se)
  }
 }
 performance[[s]]<-do.call(rbind,pp);paired[[s]]<-do.call(rbind,dd)
 screening[[s]]<-cbind(st,MeanUsedN=mean(scr$UsedN,na.rm=TRUE),MeanCIS=mean(scr$CIS,na.rm=TRUE),
  Mean_kCIS=mean(scr$kCIS,na.rm=TRUE),Mean_kRHKB=mean(scr$kRHKB,na.rm=TRUE),
  FallbackRate=mean(scr$Fallback,na.rm=TRUE),Detection=mean(scr$Detection,na.rm=TRUE),FalseExclusion=mean(scr$FalseExclusion,na.rm=TRUE))
 if(length(err))errors[[length(errors)+1L]]<-cbind(SettingID=st$SettingID,do.call(rbind,err))
 saveRDS(list(setting=st,beta_true=beta,coefficients=estimates,loss=loss,screening=scr,errors=err),
  file.path(OUTPUT_DIR,sprintf("Setting_%03d.rds",st$SettingID)))
 cat("Completed",s,"/",nrow(settings),"common successes",sum(common),"/",B,"\n")
}
comp<-do.call(rbind,comparison);perf<-do.call(rbind,performance);pair<-do.call(rbind,paired)
errtab<-if(length(errors))do.call(rbind,errors)else data.frame(SettingID=integer(),Replication=integer(),Method=character(),Message=character())
notes<-data.frame(Item=c("Methods","Design","Numerics","Interpretation","Seed ordering"),
 Detail=c("MLE; GammaRidge; RefitRidgeHKB; GammaCIS. Refit is internal only.",
 "Clean and rho=0 benchmarks retained. Same selected rows: all standardized slopes x50, response x8; no restandardization.",
 "Scaled augmented QR; MPFR reconstruction with precision escalation and agreement/residual checks; unpenalized likelihood recovery. Unresolved failures retained and logged.",
 "Common-success MSE for ranking; estimator-specific bias and variance; pairwise complete differences and pointwise Monte Carlo intervals. No guaranteed winner.",
 "Full-grid order and seed formula match the archived design. MPFR is integrated into fitting. Different R versions, tolerances and recovery policies can change exact Monte Carlo outputs; equality to the archived workbook is not asserted."))
cis_export_tables(list(Settings=settings,Comparison=comp,Performance=perf,Paired=pair,
 Screening=do.call(rbind,screening),Failures=errtab,NumericalAudit=if(length(numerical))do.call(rbind,numerical)else data.frame(Note="No shrinkage solves"),Notes=notes),file.path(OUTPUT_DIR,"Simulation_Results.xlsx"))
pdf(file.path(OUTPUT_DIR,"Simulation_Figures.pdf"),width=9,height=6)
for(sc in unique(comp$Scenario)) {
 z<-comp[comp$Scenario==sc,];v<-sapply(.cis_methods,function(m)tapply(z[[paste0("MSE_",m)]],z$rho,mean,na.rm=TRUE))
 if(all(is.finite(v))&&all(v>0)) {
  matplot(sort(unique(z$rho)),matrix(v,nrow=length(unique(z$rho))),type="b",log="y",lty=1,pch=1:4,col=.cis_colors,xlab="Predictor correlation",ylab="Mean coefficient MSE",main=sc)
  legend("topleft",.cis_methods,col=.cis_colors,lty=1,pch=1:4,bty="n")
 }
}
dev.off();writeLines(capture.output(sessionInfo()),file.path(OUTPUT_DIR,"SessionInfo.txt"))
cat("Results:",normalizePath(OUTPUT_DIR),"\n")

invisible(list(settings=settings,comparison=comp,performance=perf,paired=pair,numerical=numerical,output_dir=OUTPUT_DIR))

}
