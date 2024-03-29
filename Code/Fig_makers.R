dindex_calc<-function(data){
  dind<-pnorm(data$est,0,1)
  upper<-pnorm(data$est+1.96*data$se)
  lower<-pnorm(data$est-1.96*data$se)
  return(data.frame(Sample=data$Sample,est=dind,upper_ci=upper,lower_ci=lower))
}

dscore_mkr<-function(data,main="",colr="grey"){
  len<-dim(data)[1]
  ord=order(data$est)
  plot(data$est[ord],1:len,yaxt="n",ylab ="",xlab="",pch=21,bg=colr,xlim=c(-4,4),
       main=main,bty='n',col=colr)
  arrows(data$est[ord],1:len,data$est[ord]+1.96*data$se[ord],length = 0.05,angle = 90,col="black")
  arrows(data$est[ord],1:len,data$est[ord]-1.96*data$se[ord],length = 0.05,angle = 90,col="black")
  abline(v=0,lwd=2)
  points(data$est[ord],1:len,pch=21,bg=colr,cex=1.25)
  axis(side=2,1:len,labels = data$Sample[ord],las=T,cex.axis=.5)
  mtext("Sample (Station_Year_Rep)",side = 2,cex=1.5,padj=-6.5)
  mtext(expression("Disturbance Score " (italic("D"))),side = 1,cex=1.5,padj=3)
}

dindex_mkr<-function(data,main=""){
  len<-dim(data)[1]
  ord=order(data$est)
  plot(pnorm(data$est[ord],0,1),1:len,yaxt="n",ylab ="",xlab="",pch=21,bg="gray",xlim=c(0,1),
       main=main,bty='n')
  arrows(pnorm(data$est[ord]),1:len,pnorm(data$est[ord]+1.96*data$se[ord]),length = 0.05,angle = 90)
  arrows(pnorm(data$est[ord]),1:len,pnorm(data$est[ord]-1.96*data$se[ord]),length = 0.05,angle = 90)
  abline(v=.5,lwd=2)
  points(pnorm(data$est[ord]),1:len,pch=21,bg="gray")
  axis(side=2,1:len,labels = data$Sample[ord],las=T,cex.axis=.5)
  mtext("Station",side = 2,cex=1.5,padj=-3)
  mtext("Disturbance Index",side = 1,cex=1.5,padj=3)
}

alpha_mkr<-function(data,main="",colr="grey"){
  nsp=dim(data)[1]
  ord=order(data$est)
  par(mar=c(5,8,2,2))
  plot(data$est[ord],1:nsp,yaxt="n",ylab = "",xlab="",pch=21,bg=colr,xlim=c(-.8,.8))
  arrows(data$est[ord],1:nsp,data$est[ord]+2*data$se[ord],length = 0.05,angle = 90,col="black")
  arrows(data$est[ord],1:nsp,data$est[ord]-2*data$se[ord],length = 0.05,angle = 90,col="black")
  points(data$est[ord],1:nsp,pch=21,bg=colr,xlim=c(-.8,.8))
  
  abline(v=0,lwd=2)
  axis(side=2,1:nsp,labels = data$names[ord],las=T,cex.axis=.5)
  mtext("Taxon",side = 2,cex=1.5,padj=-7)
  mtext(bquote("Sensitivity Index ("*alpha*")"),side = 1,cex=1.5,padj=3)
}
  