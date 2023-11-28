#This code is used to make Figure 6 of the main text

#get data on taxa
df<-read.csv("Data/habitat and TOC and species abundances for Long-term 2017 and 2018 samples.csv")
sp_start<-14
sp_dat<-df[,sp_start:dim(df)[2]]
names(sp_dat)
sp_dat[sp_dat>0]<-1
#get d-scores and alphas
d_dat<-read.csv("Data/Model/D_ScoresPriors.csv")
alphas<-read.csv("Data/Model/alphas_calibration.csv")
nms<-alphas$names

#calculate richness
taxa_rich<-apply(sp_dat[,nms],1,sum)
#calculate dominance
dom<-apply((sp_dat[,nms]/apply(sp_dat[,nms],1,sum))^2,1,sum)
#calculate an ambi-like

bc_scores<-read.csv("Data/AMBI EG for 551 taxa in 2017-2018.csv")
#grab our taxa
bc_scores<-bc_scores[which(bc_scores$X551.taxa%in%nms),]
colnames(bc_scores)[1]<-"names"
#merge with alpha estimates
df_bc<-merge(alphas,bc_scores,by="names")
#remove zeros
df_bc<-df_bc[-which(df_bc$AMBI2022==0),]
wts<-c(0,1.5,3,4.5,6)
amb3<-apply(t(t(df[,df_bc$names]/apply(df[,df_bc$names],1,sum))*wts[df_bc$AMBI2022]),1,sum)

#Figure 5
pdf("Output/Comparisons.pdf")
par(mfrow=c(2,2),oma=c(2,2,2,2))
plot(df_bc$est~wts[df_bc$AMBI2022],pch=21,bg="grey",xlab="Ecological Group Weight",ylab=bquote(alpha),bty='n',cex.lab=1.2,main=paste0("corr=",round(cor(df_bc$est,wts[df_bc$AMBI2022]),3)))
mtext(text = "(a)",side=3,adj =-.3)
plot(d_dat$est~apply(sp_dat[,nms],1,sum),pch=21,bg="grey",bty="n",xlab="Taxa Richness",ylab="D-Score",cex.lab=1.2,main=paste0("corr=",round(cor(d_dat$est,taxa_rich),3)))
mtext(text = "(b)",side=3,adj = -.3)
plot(d_dat$est~amb3,pch=21,bg="grey",bty="n",xlab=bquote("AMBI"),ylab="D-Score",cex.lab=1.2,main=paste0("corr=",round(cor(d_dat$est,amb3),3)))
mtext(text = "(c)",side=3,adj = -.3)
plot(d_dat$est~dom,pch=21,bg="grey",bty="n",xlab="Simpson Index",ylab="D-Score",cex.lab=1.2,main=paste0("corr=",round(cor(d_dat$est,dom),3)))
mtext(text = "(d)",side=3,adj = -.3)
dev.off()

