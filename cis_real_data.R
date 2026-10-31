cis_bodyfat <- function(data=NULL,output_dir=file.path(getwd(),"bodyfat_Modified"),contamination_seed=123L,screen_seed=20261017L) {
for(pkg in c("MASS","openxlsx","Rmpfr"))if(!requireNamespace(pkg,quietly=TRUE))stop(paste("Install",pkg))
# -------------------------- USER SETTINGS --------------------------
OUTPUT_DIR<-output_dir
CONTAMINATION_SEED<-contamination_seed;SCREEN_SEED<-screen_seed
RATE<-.10;X_MULTIPLIER<-50;Y_MULTIPLIER<-8
if(is.null(data)) {
 if(!requireNamespace("TH.data",quietly=TRUE))stop("Supply original bodyfat data or install TH.data")
 e<-new.env();utils::data("bodyfat",package="TH.data",envir=e);raw<-e$bodyfat
} else raw<-data
stopifnot(nrow(raw)==71,all(complete.cases(raw)),all(raw$DEXfat>0))
pred<-setdiff(names(raw),"DEXfat");z<-cis_make_design(raw[,pred]);X0<-z$X;y0<-raw$DEXfat;n<-length(y0)
set.seed(CONTAMINATION_SEED);ids<-sort(sample.int(n,ceiling(RATE*n)))
X1<-X0;y1<-y0;X1[ids,-1]<-X1[ids,-1,drop=FALSE]*X_MULTIPLIER;y1[ids]<-y1[ids]*Y_MULTIPLIER
set.seed(SCREEN_SEED);before<-suppressWarnings(cis_fit_methods(X0,y0))
set.seed(SCREEN_SEED);after<-suppressWarnings(cis_fit_methods(X1,y1))
dir.create(OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
metric<-function(b) {
 if(is.null(b))return(c(GammaDeviance=NA_real_,MAE=NA_real_,RMSE=NA_real_))
 mu<-cis_positive_mean(X0,b)
 c(GammaDeviance=mean(cis_gamma_unit_deviance(y0,mu)),MAE=mean(abs(y0-mu)),RMSE=sqrt(mean((y0-mu)^2)))
}
comparison<-do.call(rbind,lapply(.cis_methods,function(m)data.frame(Method=m,t(metric(after$beta[[m]])),row.names=NULL)))
coefficients<-do.call(rbind,lapply(.cis_methods,function(m)data.frame(Method=m,Term=colnames(X0),
 Original=if(is.null(before$beta[[m]]))NA_real_ else before$beta[[m]],
 Modified=if(is.null(after$beta[[m]]))NA_real_ else after$beta[[m]])))
# Observed estimates only: no true-beta MSE or original-MLE-as-truth criterion.
vif<-function(X)diag(solve(cor(X[,-1,drop=FALSE])))
vifs<-data.frame(Variable=pred,Original=vif(X0),Modified=vif(X1))
get_screen<-function(obj,label) {
 s<-obj$meta$screen
 if(is.null(s))return(data.frame(Dataset=label,ID=seq_len(n),Injected=seq_len(n)%in%ids,RD2=NA_real_,DevianceResidual=NA_real_,Weight=NA_real_))
 data.frame(Dataset=label,ID=seq_len(n),Injected=seq_len(n)%in%ids,RD2=s$rd2,
 DevianceResidual=obj$meta$mle$rdev,Weight=s$w)
}
diagnostics<-rbind(get_screen(before,"bodyfat"),get_screen(after,"bodyfat_Modified"))
tuning<-do.call(rbind,lapply(list(Original=before,Modified=after),function(o){
 s<-o$meta$screen;data.frame(CIS=if(is.null(o$meta$CIS))NA else o$meta$CIS,
 kCIS=if(is.null(o$meta$k_CIS))NA else o$meta$k_CIS,kRHKB=if(is.null(o$meta$k_refit_HKB))NA else o$meta$k_refit_HKB,
 UsedN=if(is.null(s))NA else sum(s$w),Fallback=if(is.null(s))NA else s$fallback)}))
tuning<-cbind(Dataset=rownames(tuning),tuning);rownames(tuning)<-NULL
status<-data.frame(Method=.cis_methods,OriginalError=unname(before$errors),ModifiedError=unname(after$errors))
sc<-after$meta$screen;spectral<-data.frame();proximity<-data.frame()
if(!is.null(sc)) {
 Xr<-X1[sc$take,,drop=FALSE];ev<-eigen(crossprod(Xr*sqrt(sc$mu)),symmetric=TRUE,only.values=TRUE)$values
 spectral<-data.frame(Direction=seq_along(ev),Eigenvalue=ev,
 HKBFactor=ev/(ev+after$meta$k_refit_HKB),CISFactor=ev/(ev+after$meta$k_CIS))
}
bh<-after$beta$RefitRidgeHKB;bc<-after$beta$GammaCIS
if(!is.null(bh)&&!is.null(bc)) {
 proximity<-data.frame(Term=colnames(X0),RHKB=bh,CIS=bc,Difference=bc-bh)
}
notes<-data.frame(Item=c("Source","Perturbation","Scaling","Evaluation","Screen","Interpretation"),
 Detail=c("TH.data::bodyfat; Garcia et al. (2005), Obesity Research 13:626-634. doi:10.1038/oby.2005.67",
 paste("One fixed draw; seed",CONTAMINATION_SEED,";",length(ids),"of",n,"rows; actual percent",100*length(ids)/n),
 "Original means/SDs fixed. Same selected rows: all standardized slopes x50, y x8. Intercept stays one. No restandardization.",
 "Fit modified observations; evaluate original responses at original predictors. Restoration, not independent test prediction. No true coefficients assumed.",
 "MVE 1000 subsamples, chi-square .975; unscaled deviance cutoff 2.5. Baseline is not certified outlier-free.",
 "All three response criteria retained; coefficient differences and spectral shrinkage explain proximity without asserting a guaranteed CIS advantage."))
cis_export_tables(list(Comparison=comparison,Coefficients=coefficients,CIS_vs_RHKB=proximity,
 Spectral=spectral,VIF=vifs,Screening=diagnostics,Tuning=tuning,Status=status,
 Selected_Rows=data.frame(ID=ids),Original=raw,
 Analysis_Modified=data.frame(ID=seq_len(n),X1[,-1,drop=FALSE],DEXfat=y1),Notes=notes),
 file.path(OUTPUT_DIR,"bodyfat_Modified_Results.xlsx"))
pdf(file.path(OUTPUT_DIR,"bodyfat_Modified_Figures.pdf"),width=10,height=5)
par(mfrow=c(1,2),mar=c(5,4,3,1));barplot(vifs$Original,names.arg=pred,las=2,log="y",col=.cis_colors[4],ylab="VIF",main="Original predictor collinearity")
d<-diagnostics[diagnostics$Dataset=="bodyfat_Modified",]
if(all(is.finite(d$RD2))) {
 plot(d$ID,sqrt(d$RD2),log="y",pch=ifelse(d$Injected,17,1),col=ifelse(d$Injected,"firebrick","gray40"),xlab="Row",ylab="Robust predictor distance",main="Joint perturbation")
 abline(h=sqrt(qchisq(.975,ncol(X0)-1)),lty=2)
}
par(mfrow=c(1,3),mar=c(8,4,3,1))
for(nm in c("GammaDeviance","MAE","RMSE"))barplot(comparison[[nm]],names.arg=.cis_methods,las=2,col=.cis_colors,main=nm)
if(nrow(proximity)) {
 par(mfrow=c(1,2),mar=c(5,4,3,1));plot(proximity$RHKB,proximity$CIS,pch=19,xlab="RHKB estimated coefficients",ylab="CIS estimated coefficients",main="Coefficient proximity");abline(0,1,lty=2)
 matplot(spectral$Direction,spectral[,c("HKBFactor","CISFactor")],type="b",lty=1,pch=c(1,16),col=.cis_colors[c(3,4)],xlab="Eigen direction",ylab="Shrinkage factor",main="Realised attenuation")
 legend("bottomleft",c("RHKB","CIS"),col=.cis_colors[c(3,4)],pch=c(1,16),lty=1,bty="n")
}
dev.off();saveRDS(list(before=before,after=after,X0=X0,y0=y0,X1=X1,y1=y1,selected=ids,comparison=comparison),file.path(OUTPUT_DIR,"Fitted_Objects.rds"))
writeLines(capture.output(sessionInfo()),file.path(OUTPUT_DIR,"SessionInfo.txt"));print(comparison)
cat("Results:",normalizePath(OUTPUT_DIR),"\n")

invisible(list(comparison=comparison,coefficients=coefficients,spectral=spectral,before=before,after=after,selected=ids,output_dir=OUTPUT_DIR))

}
