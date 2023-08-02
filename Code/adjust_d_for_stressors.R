

#refit D in the while accounting for some of the potential mechanisms of D
setwd("path/to/latent_disturbance_bioassessment/")
#grab the unadjusted d-scores
d_dat<-read.csv("Data/Model/D_ScoresPriors.csv")

#grab training data
df<-read.csv("Data/habitat and TOC and species abundances for Long-term 2017 and 2018 samples.csv")
#grab sediment stressor data
sed_dat<-read.csv("Data/sediment chemistry for 2017-2018 baseline.csv")
head(sed_dat)

#grab the update stuff
scale_vars<-read.csv("Data/Model/E_calibration_means.csv")
alphas<-read.csv("Data/Model/alphas_calibration.csv")
sp_coefs<-read.csv("Data/Model/TaxaCoefficients.csv")
nms<-alphas$names


#grab species and set convert abundance to presence
sp_dat<-df[,nms]
names(sp_dat)
sp_dat[sp_dat>0]<-1

###scale test data by calibration data mean and sd
env_vars<-df[,names(df)[c(6:11,13)]]
mod_vars<-env_vars[,c("Depth","Penetration","Salinity","Temperature","Fines","Gravel")]

logit_trans<-function(x)log((x/100)/(1-x/100))

mod_vars$Fines<-logit_trans(env_vars$Fines)-(as.numeric(scale_vars[5])+as.numeric(scale_vars[6])*sqrt(env_vars$Sed_TOC))
mod_vars<-sweep(mod_vars,MARGIN = 2,STATS = as.numeric(c(scale_vars[1:4],0,scale_vars[7])),FUN = "-")
mod_vars
names(mod_vars)<-paste0(names(mod_vars),"_x")
#reassemble the scale data
df_test<-data.frame(Sample=df$Sample,sp_dat,mod_vars)

####### set up matrix for residuals
n_sites<-dim(df_test)[1]
res_new<-matrix(NA,nrow=n_sites,ncol=length(nms))
colnames(res_new)<-nms
rownames(res_new)<-df_test$Sample

#setup data for calc
i<-nms[1]
fmls<-paste0(i,'~Depth_x+Penetration_x+Salinity_x+Temperature_x+Fines_x+Gravel_x+Salinity_x:Temperature_x')
x<-model.matrix(formula(fmls),data = df_test)[,-1]

#calculate residuals
for(i in nms){
  pdct<-1/(1+exp(-cbind(rep(1,dim(df_test)[1]),x)%*%sp_coefs[,i]))
  pdct[which(pdct>9.999999e-01)]<-9.999999e-01
  a<-pbinom(df_test[,i]-1,size = 1,pdct)
  b<-dbinom(df_test[,i],size = 1,pdct)
  tmp_scr<-(2*a+b)/2
  res_new[,i]<-qlogis(tmp_scr,0,1)
}

#add the stressors to be adjusted for to the data.frame of residuals
res<-data.frame(res_new,sed_dat[match(rownames(res_new),sed_dat$Sample),c("Sed_TOC","Sed_TN")])

#set the "non-detect value of TN (i.e., 0.1) to zero
res$Sed_TN[res$Sed_TN==.1]<-0
#look a correlation of d-scores with stressors
cor(d_dat$est,res$Sed_TN)
cor(d_dat$est,res$Sed_TOC)

#set up model for lavaan calc of alphas
model_p1<-paste("D=~NA*",paste(colnames(res_new),sep=" ",collapse ="+"),"\nD~~1*D",sep="")
model<-paste(model_p1,paste(colnames(res_new),"~ Sed_TOC + Sed_TN",collapse = "\n"),"\nSed_TOC~~Sed_TN",sep = "\n")
#take a peek
cat(model)
library(lavaan)
#fit the model
fit<-sem(model,data = res,meanstructure = FALSE)
#take a peek
summary(fit)

#get the estimates
alpha_est<-fit@ParTable$est[which(fit@ParTable$op=="=~")]
#get the se 
alpha_se<-fit@ParTable$se[which(fit@ParTable$op=="=~")]

#get the TN and TOC parameters for each taxon
TN_est<-fit@ParTable$est[head(which(fit@ParTable$rhs=="Sed_TN"),-2)]
TOC_est<-fit@ParTable$est[head(which(fit@ParTable$rhs=="Sed_TOC"),-1)]
#get the se 
TN_se<-fit@ParTable$se[head(which(fit@ParTable$rhs=="Sed_TN"),-2)]
TOC_se<-fit@ParTable$se[head(which(fit@ParTable$rhs=="Sed_TOC"),-1)]
#check the we got all of them
length(TN_est)==length(nms)

#toss the info on alphas into a data.frame
alphas_new<-data.frame(names=fit@ParTable$rhs[which(fit@ParTable$op=="=~")],
                   est=alpha_est,se=alpha_se)


#set up a list for the adjusted estimates of D, D-twiddle
D_ests_new<-list()
for(i in 1:n_sites){
  ll<-function(x)sum(-log(dlogis(res_new[i,],location = alphas_new$est*x[1]+TN_est*res[i,"Sed_TN"]+TOC_est*res[i,"Sed_TOC"],scale = x[2])),na.rm = T)
  f_tem<-optim(c(0,.3),ll,hessian = T)
  D_ests_new[[i]]<-c(D=f_tem$par[1],se=sqrt(diag(solve(f_tem$hessian)))[1])
  }

#stuff these into a data frame 
d_ans<-data.frame(Sample=rownames(res),est_adj=sapply(D_ests_new,"[",1),se_adj=sapply(D_ests_new,"[",2))

#look at correlation with D before and after adjustment
cor(d_dat$est,res$Sed_TN)
cor(d_ans$est_adj,res$Sed_TN)

cor(d_dat$est,res$Sed_TOC)
cor(d_ans$est_adj,res$Sed_TOC)
#wow math works!

#grab the sediment data and merge with the estimates of d
df_stress<-merge(merge(d_dat,d_ans,by="Sample"),sed_dat,by="Sample")
names(df_stress)
stressors<-8:19
tabl<-cbind(sapply(df_stress[,stressors],function(x)round(cor.test(x,df_stress$est,use="pairwise.complete.obs")$estimate,3)),sapply(df_stress[,stressors],function(x)round(cor.test(x,df_stress$est,use="pairwise.complete.obs")$p.value,3)))
tabl[which(tabl[,2]<0.05),]

tabl2<-cbind(sapply(df_stress[,stressors],function(x)round(cor.test(x,df_stress$est_adj,use="pairwise.complete.obs")$estimate,3)),sapply(df_stress[,stressors],function(x)round(cor.test(x,df_stress$est_adj,use="pairwise.complete.obs")$p.value,3)))
tabl2[which(tabl2[,2]<0.05),]




#create smaller data.frame for plotting with shorter names for displaying the 
#adjustment of D
df_plot<-df_stress[,1:5]
nms_bits<-strsplit(df_plot$Sample,"_")
short_names<-paste(sapply(nms_bits,"[",3),sapply(nms_bits,"[",2),sapply(nms_bits,"[",4),sep = "_")

source("Code/Fig_makers.R")
df_plot$Sample<-short_names
dscore_mkr<-function(data,main="",colr="black"){
  ord<-order(data$est)
  len<-dim(data)[1]
  ord=order(data$est)
  plot(data$est[ord],1:len,yaxt="n",ylab ="",xlab="",pch=21,bg="darkgrey",xlim=c(-4,5),
       main=main,bty='n',col=colr)
  arrows(data$est[ord],1:len,data$est[ord]+1.96*data$se[ord],length = 0.1,angle = 90,col=colr)
  arrows(data$est[ord],1:len,data$est[ord]-1.96*data$se[ord],length = 0.1,angle = 90,col=colr)
  abline(v=0,lwd=2)
  points(data$est[ord],1:len,pch=21,bg=colr,cex=1.25)
  axis(side=2,1:len,labels = data$Sample[ord],las=T,cex.axis=.75)
  mtext("Sample (Station_Year_Rep)",side = 2,cex=1.5,padj=-8)
  mtext(expression("Disturbance Score" (italic("D"))),side = 1,cex=1.5,padj=3)
}

ord<-order(df_plot$est)
len<-length(ord)
pdf("Output/d_score_change_TOC_TN.pdf",width = 8,height =8)
par(oma=c(0,0,0,0),mar=c(5,10,1,1))
dscore_mkr(df_plot[,1:3],colr = "grey")
arrows(df_plot$est_adj[ord],1:len,df_plot$est_adj[ord]+1.96*df_plot$se_adj[ord],length = 0.1,angle = 90,col="black")
arrows(df_plot$est_adj[ord],1:len,df_plot$est_adj[ord]-1.96*df_plot$se_adj[ord],length = 0.1,angle = 90,col="black")
points(df_plot$est_adj[ord],1:len,pch=21,bg="white",cex=1.25)

abline(v=0,lwd=2)
legend("topleft",legend = c("D-Score",expression(tilde(D)-Score)),pch=21,pt.bg = c("grey","white"),bty='n')
dev.off()

#make plot to show the adjustment of alphas
ord_a<-order(alphas$est)
len_a<-124
pdf("Output/alpha_change_TOC_TN.pdf",width = 8,height =8)
par(oma=c(0,0,0,0),mar=c(5,9,1,1))
plot(alphas$est[ord_a],1:len_a,pch=21,bg="grey",yaxt='n',ylab="",xlab=bquote("Sensitivity Index ("*alpha*")"),xlim = c(-.8,.8),cex.lab=1.5)
arrows(alphas$est[ord_a],1:len_a,alphas$est[ord_a]+1.96*alphas$se[ord_a],length = 0.1,angle = 90,col="grey")
arrows(alphas$est[ord_a],1:len_a,alphas$est[ord_a]-1.96*alphas$se[ord_a],length = 0.1,angle = 90,col="grey")
arrows(alphas_new$est[ord_a],1:len_a,alphas_new$est[ord_a]+1.96*alphas_new$se[ord_a],length = 0.1,angle = 90,col="black")
arrows(alphas_new$est[ord_a],1:len_a,alphas_new$est[ord_a]-1.96*alphas_new$se[ord_a],length = 0.1,angle = 90,col="black")
points(alphas_new$est[ord_a],1:len_a,pch=21,bg="white")

legend("topleft",legend = c(expression(alpha),expression(tilde(alpha))),pch=21,pt.bg = c("grey","white"),bty='n')
axis(side = 2,at = 1:len_a,labels = alphas$names[ord_a],cex.axis=.5,las=T)
mtext("Taxon",side = 2,cex=1.5,padj=-7.5)
dev.off()
alphas_new
