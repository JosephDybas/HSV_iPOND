
options(scipen=999)

library(VennDiagram)
library(reshape)
library(ggplot2)
library(gplots)
library(scales)
library(RColorBrewer)
library(corrplot)
library(som)
library(cluster)
library(e1071)
library(factoextra)
##library(pcaMethods)


source("/Users/dybasj/JoeRLib/ProteomicsAnalysisFunctions.R")
source("/Users/dybasj/JoeRLib/PlottingFunctions.R")

##################
#NOT CHECKED
plotfoldchanges <- function(foldchanges1, foldchanges2, foldchanges3, figfilename) {
	pdf(figfilename)
	#combine fold change vectors
	all_foldchanges<-c(foldchanges1,foldchanges2,foldchanges3)
	#mean and standard dev
	foldchanges_mean<-mean(all_foldchanges, na.rm = TRUE)
	foldchanges_stdev<-sd(all_foldchanges, na.rm = TRUE)
	#histograms for each fold change vector
	breaksmin=floor(min(all_foldchanges,na.rm=TRUE))
	breaksmax=ceiling(max(all_foldchanges,na.rm=TRUE))
	print(breaksmin)
	print(breaksmax)
	foldchanges1_hist<-hist(foldchanges1, breaks=seq(breaksmin,breaksmax,by=0.25), plot=FALSE)
	foldchanges2_hist<-hist(foldchanges2, breaks=seq(breaksmin,breaksmax,by=0.25), plot=FALSE)
	foldchanges3_hist<-hist(foldchanges3, breaks=seq(breaksmin,breaksmax,by=0.25), plot=FALSE)
	#plot fold change histograms and add mean+/-stdev
	par(bty='l')
	plot(foldchanges1_hist$mids, foldchanges1_hist$density, type="l", col="lightgrey", lwd=3, xlab='fold changes', ylab='normalized frequency')
	par(new=T)
	plot(foldchanges2_hist$mids, foldchanges2_hist$density, type="l", col="grey", lwd=3, xlab='', ylab='', axes=FALSE)
	par(new=T)
	plot(foldchanges3_hist$mids, foldchanges3_hist$density, type="l", col="darkgrey", lwd=3, xlab='', ylab='', axes=FALSE)
	par(new=T)
	abline(v=foldchanges_mean, col="black")
	par(new=T)
	abline(v=foldchanges_mean+foldchanges_stdev, col="black")
	par(new=T)
	abline(v=foldchanges_mean-foldchanges_stdev, col="black")
	#par(new=T)
	#abline(v=1, col="salmon")
	#par(new=T)
	#abline(v=-1, col="salmon")
	dev.off()
	return
}
#
##################


##################
plotvolcanoplot <- function(normdata)
{
  
  ICP0subs<-c("P29590", "P78527", "Q16666", "P46100", "Q93009")
  #PCAsubClusIntersect<-c("Q08AF3", "P61956", "P49321", "P28340", "Q15054", "P52701", "Q6ZRQ5", "Q13111", "P29590", "P78527", "Q16666", "P46100")
  PCAsubClusIntersect<-scan("./analysisoutput_AUG2019/PCAclustering/IntersectClus_PML-PRKDC-ATRX-IFI16_20200810.txt", character(), quote = "")
  PCAsubClusUnion <-scan("./analysisoutput_AUG2019/PCAclustering/UnionClus_PML-PRKDC-ATRX-IFI16_20200810.txt", character(), quote = "")
  K11subs<-c("UBE2E1", "SMARCA2", "RPL6", "ADAR", "PLIN3", "RPL17", "RRP12")
  
  
  ##mock vs HSV WT volcano plot
	pdf(file.path(".", analysisdir, "MockHSV_VolcanoPlot_wICP0substrates_nolabels_20200810.pdf"))
	minfc <- floor(min(as.numeric(normdata$mock_hsv_avgfc), na.rm=TRUE))
	maxfc <- ceiling(max(as.numeric(normdata$mock_hsv_avgfc), na.rm=TRUE))
	minpval <- floor(min(as.numeric(normdata$mock_hsv_neglog10pval), na.rm=TRUE))
	maxpval <- ceiling(max(as.numeric(normdata$mock_hsv_neglog10pval), na.rm=TRUE))
	#if ( abs(minfc)>maxfc ){
	#	maxfc<-abs(minfc)
	#}else {
	#	minfc <- -maxfc
	#}
	#minfc <- -6
	#maxfc <- 5
	#minpval <- 0
	#maxpval <- 3.5
	par(bty='l',lwd=3)
	plot(normdata$mock_hsv_avgfc, normdata$mock_hsv_neglog10pval,    
	     pch=20,
	     cex=1,
	     #col="grey50",
	     col=ifelse(normdata$mock_hsv_neglog10pval<1.3, "grey75", "grey50"),
	     xlim=c(minfc,maxfc), ylim=c(minpval,maxpval),
	     xlab="",ylab="",xaxt='n',yaxt='n'
  )
	#
	#add known substrates
	for ( i in 1:length(ICP0subs) ){
	  currsub<-ICP0subs[i]
	#  #points(normdata[currsub,"mock_hsv_avgfc"],normdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.25, col="red")
	  text(normdata[currsub,"mock_hsv_avgfc"],normdata[currsub,"mock_hsv_neglog10pval"], labels=normdata[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
	}
	##
	#add clustered proteins
	for ( i in 1:length(PCAsubClusUnion) ){
	  currsub<-PCAsubClusUnion[i]
	  points(normdata[currsub,"mock_hsv_avgfc"],normdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.25, col="dodgerblue")
	#  text(normdata[currsub,"mock_hsv_avgfc"],normdata[currsub,"mock_hsv_neglog10pval"], labels=normdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
	}
	#
	#add clustered proteins
	#for ( i in 1:length(PCAsubClusIntersect) ){
	#  currsub<-PCAsubClusIntersect[i]
	#  points(normdata[currsub,"mock_hsv_avgfc"],normdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.25, col="navy")
	#  #text(normdata[currsub,"mock_hsv_avgfc"],normdata[currsub,"mock_hsv_neglog10pval"], labels=normdata[currsub,"geneID"], cex=0.75, pos=2, col="navy", font=2)
	#}
	#add K11 subs
	#for ( i in 1:length(K11subs) ){
	#  currsub<-K11subs[i]
	#  if ( currsub %in% normdata$geneID  ){
	#    points(normdata[normdata$geneID==currsub,"mock_hsv_avgfc"],normdata[normdata$geneID==currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.25, col="olivedrab")
	#    text(normdata[normdata$geneID==currsub,"mock_hsv_avgfc"],normdata[normdata$geneID==currsub,"mock_hsv_neglog10pval"], labels=normdata[normdata$geneID==currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
	#  }
	#}
	#
	axis(1,c(seq(minfc,maxfc,by=1)),col="black",cex.axis=1,lwd=3)
	mtext("HSV-mock log2 fold change",1,line=2.5)
	axis(2,c(seq(minpval,maxpval,by=0.5)),col="black",cex.axis=1,lwd=3)
	mtext("HSV-mock -log10 pval",2,line=2.5)
	dev.off()

	##mock vs dICP0 volcano plot
	pdf(file.path(".", analysisdir, "MockdICP0_VolcanoPlot_wICP0substrates_nolabels_20200810.pdf"))
	minfc <- floor(min(as.numeric(normdata$mock_dICP0_avgfc), na.rm=TRUE))
	maxfc <- ceiling(max(as.numeric(normdata$mock_dICP0_avgfc), na.rm=TRUE))
	minpval <- floor(min(as.numeric(normdata$mock_dICP0_neglog10pval), na.rm=TRUE))
	maxpval <- ceiling(max(as.numeric(normdata$mock_dICP0_neglog10pval), na.rm=TRUE))
	#if ( abs(minfc)>maxfc ){
	#	maxfc<-abs(minfc)
	#}else {
	#	minfc<--maxfc
	#}
	#minfc <- -5
	#maxfc <- 5
	#minpval <- 0
	#maxpval <- 3
	par(bty='l',lwd=3)
	plot(normdata$mock_dICP0_avgfc, normdata$mock_dICP0_neglog10pval,
	     pch=20,
	     cex=1,
	     #col="grey50",
	     col=ifelse(normdata$mock_dICP0_neglog10pval<1.3, "grey75", "grey50"),
	     xlim=c(minfc,maxfc), ylim=c(minpval,maxpval),
	     xlab="",ylab="",xaxt='n',yaxt='n'
	)
	#
	#add known substrates
	for ( i in 1:length(ICP0subs) ){
	  currsum<-ICP0subs[i]
	#  points(normdata[currsum,"mock_dICP0_avgfc"],normdata[currsum,"mock_dICP0_neglog10pval"], pch=20, cex=1.25, col="red")
	  text(normdata[currsum,"mock_dICP0_avgfc"],normdata[currsum,"mock_dICP0_neglog10pval"], labels=normdata[currsum,"geneID"], cex=0.75, pos=2, col="red", font=2)
	}
	#
	#add clustered proteins
	for ( i in 1:length(PCAsubClusUnion) ){
	  currsub<-PCAsubClusUnion[i]
	  points(normdata[currsub,"mock_dICP0_avgfc"],normdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.25, col="dodgerblue")
	#  text(normdata[currsub,"mock_dICP0_avgfc"],normdata[currsub,"mock_dICP0_neglog10pval"], labels=normdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
	}
	#
	#add clustered proteins
	#for ( i in 1:length(PCAsubClusIntersect) ){
	#  currsub<-PCAsubClusIntersect[i]
	#  points(normdata[currsub,"mock_dICP0_avgfc"],normdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.25, col="navy")
	#  text(normdata[currsub,"mock_dICP0_avgfc"],normdata[currsub,"mock_dICP0_neglog10pval"], labels=normdata[currsub,"geneID"], cex=0.75, pos=2, col="navy", font=2)
	#}
	#add K11 subs
	#for ( i in 1:length(K11subs) ){
	#  currsub<-K11subs[i]
	#  if ( currsub %in% normdata$geneID  ){
	#    points(normdata[normdata$geneID==currsub,"mock_dICP0_avgfc"],normdata[normdata$geneID==currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.25, col="olivedrab")
	#    text(normdata[normdata$geneID==currsub,"mock_dICP0_avgfc"],normdata[normdata$geneID==currsub,"mock_dICP0_neglog10pval"], labels=normdata[normdata$geneID==currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
	#  }
	#}
	#
	axis(1,c(seq(minfc,maxfc,by=1)),col="black",cex.axis=1,lwd=3)
	mtext("dICP0-mock log2 fold change",1,line=2.5)
	axis(2,c(seq(minpval,maxpval,by=0.5)),col="black",cex.axis=1,lwd=3)
	mtext("dICP0-mock -log10 pval",2,line=2.5)
	dev.off()
	
	##HSV WT vs dICP0 volcano plot
	pdf(file.path(".", analysisdir, "HSVdICP0_VolcanoPlot_wICP0substrates_nolabels_20200810.pdf"))
	minfc <- floor(min(as.numeric(normdata$hsv_dICP0_avgfc), na.rm=TRUE))
	maxfc <- ceiling(max(as.numeric(normdata$hsv_dICP0_avgfc), na.rm=TRUE))
	minpval <- floor(min(as.numeric(normdata$hsv_dICP0_neglog10pval), na.rm=TRUE))
	maxpval <- ceiling(max(as.numeric(normdata$hsv_dICP0_neglog10pval), na.rm=TRUE))
	#if ( abs(minfc)>maxfc ){
	#	maxfc <- abs(minfc)
	#}else {
	#	minfc <- -maxfc
	#}
	#minfc <- -6
	#maxfc <- 6
	#minpval <- 0
	#maxpval <- 3
	par(bty='l',lwd=3)
	plot(normdata$hsv_dICP0_avgfc, normdata$hsv_dICP0_neglog10pval,
	     pch=20,
	     cex=1,
	     #col="grey50",
	     col=ifelse(normdata$hsv_dICP0_neglog10pval<1.3, "grey75", "grey50"),
	     xlim=c(minfc,maxfc), ylim=c(minpval,maxpval),
	     xlab="",ylab="",xaxt='n',yaxt='n'
	)
	#
	#add known substrates
	for ( i in 1:length(ICP0subs) ){
	  currsum<-ICP0subs[i]
	#  points(normdata[currsum,"hsv_dICP0_avgfc"],normdata[currsum,"hsv_dICP0_neglog10pval"], pch=20, cex=1.25, col="red")
	  text(normdata[currsum,"hsv_dICP0_avgfc"],normdata[currsum,"hsv_dICP0_neglog10pval"], labels=normdata[currsum,"geneID"], cex=0.75, pos=2, col="red", font=2)
	}
	#
	#add clustered proteins
	for ( i in 1:length(PCAsubClusUnion) ){
	  currsub<-PCAsubClusUnion[i]
	  points(normdata[currsub,"hsv_dICP0_avgfc"],normdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.25, col="dodgerblue")
	#  text(normdata[currsub,"hsv_dICP0_avgfc"],normdata[currsub,"hsv_dICP0_neglog10pval"], labels=normdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
	}
	#
	#add clustered proteins
	#for ( i in 1:length(PCAsubClusIntersect) ){
	#  currsub<-PCAsubClusIntersect[i]
	#  points(normdata[currsub,"hsv_dICP0_avgfc"],normdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.25, col="navy")
	#  text(normdata[currsub,"hsv_dICP0_avgfc"],normdata[currsub,"hsv_dICP0_neglog10pval"], labels=normdata[currsub,"geneID"], cex=0.75, pos=2, col="navy", font=2)
	#}
	#add K11 subs
	#for ( i in 1:length(K11subs) ){
	#  currsub<-K11subs[i]
	#  if ( currsub %in% normdata$geneID  ){
	#    points(normdata[normdata$geneID==currsub,"hsv_dICP0_avgfc"],normdata[normdata$geneID==currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.25, col="olivedrab")
	#    text(normdata[normdata$geneID==currsub,"hsv_dICP0_avgfc"],normdata[normdata$geneID==currsub,"hsv_dICP0_neglog10pval"], labels=normdata[normdata$geneID==currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
	#  }
	#}
	#
	axis(1,c(seq(minfc,maxfc,by=1)),col="black",cex.axis=1,lwd=3)
	mtext("dICP0-HSV log2 fold change",1,line=2.5)
	axis(2,c(seq(minpval,maxpval,by=0.5)),col="black",cex.axis=1,lwd=3)
	mtext("dICP0-HSV -log10 pval",2,line=2.5)
	dev.off()
	
	return()

}
#
##################


if ( 0 )
{
plotfcabundance <- function(repfc,normdata_zscores,expsigtest)
{

	#####clusteredproteins <- scan("iPONDdataPCA_Abundance_PCA123ClusterUnion_TopFoldChanges.txt", what="")
	#####clusteredproteins_num<-length(clusteredproteins)
	
	#pdf("MockHSV_FoldChangeAbundancePlot_wICP0substrates.pdf")
	pdf("MockHSV_FoldChangeAbundancePlot_wICP0substrateClusters_20171017.pdf")
	minfc<-floor(min(repfc$mock_hsv_avgfc))
	maxfc<-ceiling(max(repfc$mock_hsv_avgfc))
	if ( abs(minfc)>maxfc ){
		maxfc<-abs(minfc)
	}else {
		minfc<--maxfc
	}
	minz<-floor(min(normdata_zscores$hsv_avgzscore))
	maxz<-ceiling(max(normdata_zscores$hsv_avgzscore))
	par(bty='l')
	proteins<-intersect(row.names(repfc),row.names(normdata_zscores))
	plot(repfc[proteins,"mock_hsv_avgfc"], normdata_zscores[proteins,"hsv_avgzscore"],
	     pch=20,
		 cex=0.75,
		 col=ifelse(expsigtest[proteins,"mock_hsv_neglog10pval"]<1, "lightgrey", ifelse(abs(repfc[proteins,"mock_hsv_avgfc"])<1, "lightskyblue", ifelse(abs(repfc[proteins,"mock_hsv_avgfc"])>1.585, "royalblue4", "royalblue1"))),
		 xlim=c(-8,8), ylim=c(-3,5),
		 xlab="",ylab="",xaxt='n',yaxt='n')
	addproteins<-row.names(subset(expsigtest, expsigtest[proteins,"mock_hsv_neglog10pval"]>1))
	points(repfc[addproteins,"mock_hsv_avgfc"], normdata_zscores[addproteins,"hsv_avgzscore"],
	       pch=20,
		   cex=0.75,
		   col=ifelse(abs(repfc[addproteins,"mock_hsv_avgfc"])<1, "lightskyblue", ifelse(abs(repfc[addproteins,"mock_hsv_avgfc"])>1.585, "royalblue4", "royalblue1")))
	abline(h=0,v=0)
	
	######add proteins clustered with known substrates
	#####for ( clusprot in clusteredproteins ){
	#####	points(repfc[clusprot,"mock_hsv_avgfc"],normdata_zscores[clusprot,"hsv_avgzscore"], pch=20, col="salmon")
	#####}

	#add known substrates
	points(repfc["P29590","mock_hsv_avgfc"],normdata_zscores["P29590","hsv_avgzscore"], pch=20, col="red")
	text(repfc["P29590","mock_hsv_avgfc"],normdata_zscores["P29590","hsv_avgzscore"], labels="PML", cex=0.75, pos=2, col="red", font=2)
	points(repfc["P78527","mock_hsv_avgfc"],normdata_zscores["P78527","hsv_avgzscore"], pch=20, col="red")
	text(repfc["P78527","mock_hsv_avgfc"],normdata_zscores["P78527","hsv_avgzscore"], labels="PRKDC", cex=0.75, pos=2, col="red", font=2)
	points(repfc["Q16666","mock_hsv_avgfc"],normdata_zscores["Q16666","hsv_avgzscore"], pch=20, col="red")
	text(repfc["Q16666","mock_hsv_avgfc"],normdata_zscores["Q16666","hsv_avgzscore"], labels="IFI16", cex=0.75, pos=3, col="red", font=2)
	points(repfc["Q93009","mock_hsv_avgfc"],normdata_zscores["Q93009","hsv_avgzscore"], pch=20, col="red")
	text(repfc["Q93009","mock_hsv_avgfc"],normdata_zscores["Q93009","hsv_avgzscore"], labels="USP7", cex=0.75, pos=3, col="red", font=2)
	
	#SLFN5 Q08AF3
	points(repfc["Q08AF3","mock_hsv_avgfc"],normdata_zscores["Q08AF3","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q08AF3","mock_hsv_avgfc"],normdata_zscores["Q08AF3","hsv_avgzscore"], labels="SLFN5", cex=0.75, pos=2, col="red", font=2)
	#SUMO2 P61956
	points(repfc["P61956","mock_hsv_avgfc"],normdata_zscores["P61956","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P61956","mock_hsv_avgfc"],normdata_zscores["P61956","hsv_avgzscore"], labels="SUMO2", cex=0.75, pos=2, col="red", font=2)
	#NASP P49321
	points(repfc["P49321","mock_hsv_avgfc"],normdata_zscores["P49321","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P49321","mock_hsv_avgfc"],normdata_zscores["P49321","hsv_avgzscore"], labels="NASP", cex=0.75, pos=2, col="red", font=2)
	#POLD1 P28340
	points(repfc["P28340","mock_hsv_avgfc"],normdata_zscores["P28340","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P28340","mock_hsv_avgfc"],normdata_zscores["P28340","hsv_avgzscore"], labels="POLD1", cex=0.75, pos=2, col="red", font=2)
	#POLD3 Q15054
	points(repfc["Q15054","mock_hsv_avgfc"],normdata_zscores["Q15054","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q15054","mock_hsv_avgfc"],normdata_zscores["Q15054","hsv_avgzscore"], labels="POLD3", cex=0.75, pos=2, col="red", font=2)
	#MSH6 P52701
	points(repfc["P52701","mock_hsv_avgfc"],normdata_zscores["P52701","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P52701","mock_hsv_avgfc"],normdata_zscores["P52701","hsv_avgzscore"], labels="MSH6", cex=0.75, pos=2, col="red", font=2)
	#MMS22L Q6ZRQ5
	points(repfc["Q6ZRQ5","mock_hsv_avgfc"],normdata_zscores["Q6ZRQ5","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q6ZRQ5","mock_hsv_avgfc"],normdata_zscores["Q6ZRQ5","hsv_avgzscore"], labels="MMS22L", cex=0.75, pos=2, col="red", font=2)
	#CHAF1A Q13111
	points(repfc["Q13111","mock_hsv_avgfc"],normdata_zscores["Q13111","hsv_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q13111","mock_hsv_avgfc"],normdata_zscores["Q13111","hsv_avgzscore"], labels="CHAF1A", cex=0.75, pos=2, col="red", font=2)
	
	#substrate predictions based on zscore and fold change patterns
	#points(repfc["Q13310","mock_hsv_avgfc"],normdata_zscores["Q13310","hsv_avgzscore"], pch=20, col="red")
	#text(repfc["Q13310","mock_hsv_avgfc"],normdata_zscores["Q13310","hsv_avgzscore"], labels="PABPC4", cex=0.75, pos=1, col="red", font=2)
	#points(repfc["Q5SSJ5","mock_hsv_avgfc"],normdata_zscores["Q5SSJ5","hsv_avgzscore"], pch=20, col="red")
	#text(repfc["Q5SSJ5","mock_hsv_avgfc"],normdata_zscores["Q5SSJ5","hsv_avgzscore"], labels="HP1BP3", cex=0.75, pos=1, col="red", font=2)
	#points(repfc["O15371","mock_hsv_avgfc"],normdata_zscores["O15371","hsv_avgzscore"], pch=20, col="red")
	#text(repfc["O15371","mock_hsv_avgfc"],normdata_zscores["O15371","hsv_avgzscore"], labels="EIF3D", cex=0.75, pos=2, col="red", font=2)
	#points(repfc["P17600","mock_hsv_avgfc"],normdata_zscores["P17600","hsv_avgzscore"], pch=20, col="red")
	#text(repfc["P17600","mock_hsv_avgfc"],normdata_zscores["P17600","hsv_avgzscore"], labels="SYN1", cex=0.75, pos=3, col="red", font=2)
	#points(repfc["Q9UQ80","mock_hsv_avgfc"],normdata_zscores["Q9UQ80","hsv_avgzscore"], pch=20, col="red")
	#text(repfc["Q9UQ80","mock_hsv_avgfc"],normdata_zscores["Q9UQ80","hsv_avgzscore"], labels="PA2G4", cex=0.75, pos=1, col="red", font=2)
	
	axis(1,c(seq(-8,8,by=1)),col="black",cex.axis=0.75)
	mtext("HSV-mock log2 fold change",1,line=2.5)
	axis(2,c(seq(-3,5,by=1)),col="black",cex.axis=0.75)
	mtext("HSV abundance z-score",2,line=2.5)
	dev.off()

	#pdf("MockdICP0_FoldChangeAbundancePlot_wICP0substrates.pdf")
	pdf("MockdICP0_FoldChangeAbundancePlot_wICP0substrateClusters_20171017.pdf")
	minfc<-floor(min(repfc$mock_dICP0_avgfc))
	maxfc<-ceiling(max(repfc$mock_dICP0_avgfc))
	if ( abs(minfc)>maxfc ){
		maxfc<-abs(minfc)
	}else {
		minfc<--maxfc
	}
	minz<-floor(min(normdata_zscores$dICP0_avgzscore))
	maxz<-ceiling(max(normdata_zscores$dICP0_avgzscore))
	par(bty='l')
	proteins<-intersect(row.names(repfc),row.names(normdata_zscores))
	plot(repfc[proteins,"mock_dICP0_avgfc"], normdata_zscores[proteins,"dICP0_avgzscore"],
         pch=20,
	     cex=0.75,
		 col=ifelse(expsigtest[proteins,"mock_dICP0_neglog10pval"]<1, "lightgrey", ifelse(abs(repfc[proteins,"mock_dICP0_avgfc"])<1, "lightskyblue", ifelse(abs(repfc[proteins,"mock_dICP0_avgfc"])>1.585, "royalblue4", "royalblue1"))),
		 xlim=c(-8,8), ylim=c(-3,5),
	     xlab="",ylab="",xaxt='n',yaxt='n')
	addproteins<-row.names(subset(expsigtest, expsigtest[proteins,"mock_dICP0_neglog10pval"]>1))
	points(repfc[addproteins,"mock_dICP0_avgfc"], normdata_zscores[addproteins,"dICP0_avgzscore"],
	       pch=20,
		   cex=0.75,
		   col=ifelse(abs(repfc[addproteins,"mock_dICP0_avgfc"])<1, "lightskyblue", ifelse(abs(repfc[addproteins,"mock_dICP0_avgfc"])>1.585, "royalblue4", "royalblue1")))
	abline(h=0,v=0)
	
	######add proteins clustered with known substrates
	#####for ( clusprot in clusteredproteins ){
	#####	points(repfc[clusprot,"mock_dICP0_avgfc"],normdata_zscores[clusprot,"dICP0_avgzscore"], pch=20, col="salmon")
	#####}
	
	#add known substrates
	points(repfc["P29590","mock_dICP0_avgfc"],normdata_zscores["P29590","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["P29590","mock_dICP0_avgfc"],normdata_zscores["P29590","dICP0_avgzscore"], labels="PML", cex=0.75, pos=2, col="red", font=2)
	points(repfc["P78527","mock_dICP0_avgfc"],normdata_zscores["P78527","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["P78527","mock_dICP0_avgfc"],normdata_zscores["P78527","dICP0_avgzscore"], labels="PRKDC", cex=0.75, pos=2, col="red", font=2)
	points(repfc["Q16666","mock_dICP0_avgfc"],normdata_zscores["Q16666","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["Q16666","mock_dICP0_avgfc"],normdata_zscores["Q16666","dICP0_avgzscore"], labels="IFI16", cex=0.75, pos=2, col="red", font=2)
	points(repfc["Q93009","mock_dICP0_avgfc"],normdata_zscores["Q93009","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["Q93009","mock_dICP0_avgfc"],normdata_zscores["Q93009","dICP0_avgzscore"], labels="USP7", cex=0.75, pos=3, col="red", font=2)
	
	#SLFN5 Q08AF3
	points(repfc["Q08AF3","mock_dICP0_avgfc"],normdata_zscores["Q08AF3","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q08AF3","mock_dICP0_avgfc"],normdata_zscores["Q08AF3","dICP0_avgzscore"], labels="SLFN5", cex=0.75, pos=2, col="red", font=2)
	#SUMO2 P61956
	points(repfc["P61956","mock_dICP0_avgfc"],normdata_zscores["P61956","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P61956","mock_dICP0_avgfc"],normdata_zscores["P61956","dICP0_avgzscore"], labels="SUMO2", cex=0.75, pos=2, col="red", font=2)
	#NASP P49321
	points(repfc["P49321","mock_dICP0_avgfc"],normdata_zscores["P49321","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P49321","mock_dICP0_avgfc"],normdata_zscores["P49321","dICP0_avgzscore"], labels="NASP", cex=0.75, pos=2, col="red", font=2)
	#POLD1 P28340
	points(repfc["P28340","mock_dICP0_avgfc"],normdata_zscores["P28340","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P28340","mock_dICP0_avgfc"],normdata_zscores["P28340","dICP0_avgzscore"], labels="POLD1", cex=0.75, pos=2, col="red", font=2)
	#POLD3 Q15054
	points(repfc["Q15054","mock_dICP0_avgfc"],normdata_zscores["Q15054","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q15054","mock_dICP0_avgfc"],normdata_zscores["Q15054","dICP0_avgzscore"], labels="POLD3", cex=0.75, pos=2, col="red", font=2)
	#MSH6 P52701
	points(repfc["P52701","mock_dICP0_avgfc"],normdata_zscores["P52701","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P52701","mock_dICP0_avgfc"],normdata_zscores["P52701","dICP0_avgzscore"], labels="MSH6", cex=0.75, pos=2, col="red", font=2)
	#MMS22L Q6ZRQ5
	points(repfc["Q6ZRQ5","mock_dICP0_avgfc"],normdata_zscores["Q6ZRQ5","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q6ZRQ5","mock_dICP0_avgfc"],normdata_zscores["Q6ZRQ5","dICP0_avgzscore"], labels="MMS22L", cex=0.75, pos=2, col="red", font=2)
	#CHAF1A Q13111
	points(repfc["Q13111","mock_dICP0_avgfc"],normdata_zscores["Q13111","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q13111","mock_dICP0_avgfc"],normdata_zscores["Q13111","dICP0_avgzscore"], labels="CHAF1A", cex=0.75, pos=2, col="red", font=2)
	
	#substrate predictions based on zscore and fold change patterns
	#points(repfc["Q13310","mock_dICP0_avgfc"],normdata_zscores["Q13310","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["Q13310","mock_dICP0_avgfc"],normdata_zscores["Q13310","dICP0_avgzscore"], labels="PABPC4", cex=0.75, pos=1, col="red", font=2)
	#points(repfc["Q5SSJ5","mock_dICP0_avgfc"],normdata_zscores["Q5SSJ5","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["Q5SSJ5","mock_dICP0_avgfc"],normdata_zscores["Q5SSJ5","dICP0_avgzscore"], labels="HP1BP3", cex=0.75, pos=3, col="red", font=2)
	#points(repfc["O15371","mock_dICP0_avgfc"],normdata_zscores["O15371","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["O15371","mock_dICP0_avgfc"],normdata_zscores["O15371","dICP0_avgzscore"], labels="EIF3D", cex=0.75, pos=2, col="red", font=2)
	#points(repfc["P17600","mock_dICP0_avgfc"],normdata_zscores["P17600","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["P17600","mock_dICP0_avgfc"],normdata_zscores["P17600","dICP0_avgzscore"], labels="SYN1", cex=0.75, pos=4, col="red", font=2)
	#points(repfc["Q9UQ80","mock_dICP0_avgfc"],normdata_zscores["Q9UQ80","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["Q9UQ80","mock_dICP0_avgfc"],normdata_zscores["Q9UQ80","dICP0_avgzscore"], labels="PA2G4", cex=0.75, pos=4, col="red", font=2)
	
	axis(1,c(seq(-8,8,by=1)),col="black",cex.axis=0.75)
	mtext("dIPC0-mock log2 fold change",1,line=2.5)
	axis(2,c(seq(-3,5,by=1)),col="black",cex.axis=0.75)
	mtext("dICP0 abundance z-score",2,line=2.5)
	dev.off()

	#pdf("HSVdICP0_FoldChangeAbundancePlot_wICP0substrates.pdf")
	pdf("HSVdICP0_FoldChangeAbundancePlot_wICP0substrateClusters_20171017.pdf")
	minfc<-floor(min(repfc$hsv_dICP0_avgfc))
	maxfc<-ceiling(max(repfc$hsv_dICP0_avgfc))
	if ( abs(minfc)>maxfc ){
		maxfc<-abs(minfc)
	}else {
		minfc<--maxfc
	}
	minz<-floor(min(normdata_zscores$dICP0_avgzscore))
	maxz<-ceiling(max(normdata_zscores$dICP0_avgzscore))
	par(bty='l')
	proteins<-intersect(row.names(repfc),row.names(normdata_zscores))
	plot(repfc[proteins,"hsv_dICP0_avgfc"], normdata_zscores[proteins,"dICP0_avgzscore"],
         pch=20,
	     cex=0.75,
		 col=ifelse(expsigtest[proteins,"hsv_dICP0_neglog10pval"]<1, "lightgrey", ifelse(abs(repfc[proteins,"hsv_dICP0_avgfc"])<1, "lightskyblue", ifelse(abs(repfc[proteins,"hsv_dICP0_avgfc"])>1.585, "royalblue4", "royalblue1"))),
		 xlim=c(-8,8), ylim=c(-3,5),
	     xlab="",ylab="",xaxt='n',yaxt='n')
	addproteins<-row.names(subset(expsigtest, expsigtest[proteins,"hsv_dICP0_neglog10pval"]>1))
	points(repfc[addproteins,"hsv_dICP0_avgfc"], normdata_zscores[addproteins,"dICP0_avgzscore"],
	       pch=20,
		   cex=0.75,
		   col=ifelse(abs(repfc[addproteins,"hsv_dICP0_avgfc"])<1, "lightskyblue", ifelse(abs(repfc[addproteins,"hsv_dICP0_avgfc"])>1.585, "royalblue4", "royalblue1")))
	abline(h=0,v=0)
	
	######add proteins clustered with known substrates
	#####for ( clusprot in clusteredproteins ){
	#####	points(repfc[clusprot,"hsv_dICP0_avgfc"],normdata_zscores[clusprot,"dICP0_avgzscore"], pch=20, col="salmon")
	#####}
	
	#add known substrates
	points(repfc["P29590","hsv_dICP0_avgfc"],normdata_zscores["P29590","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["P29590","hsv_dICP0_avgfc"],normdata_zscores["P29590","dICP0_avgzscore"], labels="PML", cex=0.75, pos=2, col="red", font=2)
	points(repfc["P78527","hsv_dICP0_avgfc"],normdata_zscores["P78527","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["P78527","hsv_dICP0_avgfc"],normdata_zscores["P78527","dICP0_avgzscore"], labels="PRKDC", cex=0.75, pos=2, col="red", font=2)
	points(repfc["Q16666","hsv_dICP0_avgfc"],normdata_zscores["Q16666","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["Q16666","hsv_dICP0_avgfc"],normdata_zscores["Q16666","dICP0_avgzscore"], labels="IFI16", cex=0.75, pos=2, col="red", font=2)
	points(repfc["Q93009","hsv_dICP0_avgfc"],normdata_zscores["Q93009","dICP0_avgzscore"], pch=20, col="red")
	text(repfc["Q93009","hsv_dICP0_avgfc"],normdata_zscores["Q93009","dICP0_avgzscore"], labels="USP7", cex=0.75, pos=4, col="red", font=2)
	
	#SLFN5 Q08AF3
	points(repfc["Q08AF3","hsv_dICP0_avgfc"],normdata_zscores["Q08AF3","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q08AF3","hsv_dICP0_avgfc"],normdata_zscores["Q08AF3","dICP0_avgzscore"], labels="SLFN5", cex=0.75, pos=2, col="red", font=2)
	#SUMO2 P61956
	points(repfc["P61956","hsv_dICP0_avgfc"],normdata_zscores["P61956","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P61956","hsv_dICP0_avgfc"],normdata_zscores["P61956","dICP0_avgzscore"], labels="SUMO2", cex=0.75, pos=2, col="red", font=2)
	#NASP P49321
	points(repfc["P49321","hsv_dICP0_avgfc"],normdata_zscores["P49321","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P49321","hsv_dICP0_avgfc"],normdata_zscores["P49321","dICP0_avgzscore"], labels="NASP", cex=0.75, pos=2, col="red", font=2)
	#POLD1 P28340
	points(repfc["P28340","hsv_dICP0_avgfc"],normdata_zscores["P28340","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P28340","hsv_dICP0_avgfc"],normdata_zscores["P28340","dICP0_avgzscore"], labels="POLD1", cex=0.75, pos=2, col="red", font=2)
	#POLD3 Q15054
	points(repfc["Q15054","hsv_dICP0_avgfc"],normdata_zscores["Q15054","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q15054","hsv_dICP0_avgfc"],normdata_zscores["Q15054","dICP0_avgzscore"], labels="POLD3", cex=0.75, pos=2, col="red", font=2)
	#MSH6 P52701
	points(repfc["P52701","hsv_dICP0_avgfc"],normdata_zscores["P52701","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["P52701","hsv_dICP0_avgfc"],normdata_zscores["P52701","dICP0_avgzscore"], labels="MSH6", cex=0.75, pos=2, col="red", font=2)
	#MMS22L Q6ZRQ5
	points(repfc["Q6ZRQ5","hsv_dICP0_avgfc"],normdata_zscores["Q6ZRQ5","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q6ZRQ5","hsv_dICP0_avgfc"],normdata_zscores["Q6ZRQ5","dICP0_avgzscore"], labels="MMS22L", cex=0.75, pos=2, col="red", font=2)
	#CHAF1A Q13111
	points(repfc["Q13111","hsv_dICP0_avgfc"],normdata_zscores["Q13111","dICP0_avgzscore"], pch=20, cex=2, col="red")
	text(repfc["Q13111","hsv_dICP0_avgfc"],normdata_zscores["Q13111","dICP0_avgzscore"], labels="CHAF1A", cex=0.75, pos=2, col="red", font=2)
	
	#substrate predictions based on zscore and fold change patterns
	#points(repfc["Q13310","hsv_dICP0_avgfc"],normdata_zscores["Q13310","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["Q13310","hsv_dICP0_avgfc"],normdata_zscores["Q13310","dICP0_avgzscore"], labels="PABPC4", cex=0.75, pos=1, col="red", font=2)
	#points(repfc["Q5SSJ5","hsv_dICP0_avgfc"],normdata_zscores["Q5SSJ5","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["Q5SSJ5","hsv_dICP0_avgfc"],normdata_zscores["Q5SSJ5","dICP0_avgzscore"], labels="HP1BP3", cex=0.75, pos=4, col="red", font=2)
	#points(repfc["O15371","hsv_dICP0_avgfc"],normdata_zscores["O15371","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["O15371","hsv_dICP0_avgfc"],normdata_zscores["O15371","dICP0_avgzscore"], labels="EIF3D", cex=0.75, pos=2, col="red", font=2)
	#points(repfc["P17600","hsv_dICP0_avgfc"],normdata_zscores["P17600","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["P17600","hsv_dICP0_avgfc"],normdata_zscores["P17600","dICP0_avgzscore"], labels="SYN1", cex=0.75, pos=3, col="red", font=2)
	#points(repfc["Q9UQ80","hsv_dICP0_avgfc"],normdata_zscores["Q9UQ80","dICP0_avgzscore"], pch=20, col="red")
	#text(repfc["Q9UQ80","hsv_dICP0_avgfc"],normdata_zscores["Q9UQ80","dICP0_avgzscore"], labels="PA2G4", cex=0.75, pos=3, col="red", font=2)
	
	axis(1,c(seq(-8,8,by=1)),col="black",cex.axis=0.75)
	mtext("dIPC0-HSV log2 fold change",1,line=2.5)
	axis(2,c(seq(-3,5,by=1)),col="black",cex.axis=0.75)
	mtext("dICP0 abundance z-score",2,line=2.5)
	dev.off()
	
	return()

}
}

snakeplot <- function(normdata)
{
  
  pdf("Mock_SnakePlot_wICP0substrates_20171017.pdf")
  mockzscoredata<-normdata[,"mock_avgzscore"]
  normdata_sorted<-normdata[order(normdata$mock_avgzscore, decreasing=TRUE), ]
  normdata_sorted$rownum<-1:nrow(normdata_sorted)
  minzscore<-floor(min(normdata_sorted$mock_avgzscore, na.rm=TRUE))
  maxzscore<-ceiling(max(normdata_sorted$mock_avgzscore, na.rm=TRUE))
  #if ( abs(minzscore)>maxzscore ){
  #	maxzscore<-abs(minzscore)
  #}else {
  #	minzscore <- -maxzscore
  #}
	par(bty='l',lwd=3)
	plot(normdata_sorted$mock_avgzscore,
	     pch=20,
	     cex=1.5,
	     col="grey70",
	     ylim=c(minzscore,maxzscore),
	     xlab="", ylab="", xaxt='n', yaxt='n')
	
	#PML
	points(normdata_sorted["P29590","rownum"], normdata_sorted["P29590","mock_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P29590","rownum"], normdata_sorted["P29590","mock_avgzscore"], labels="PML", cex=0.75, pos=4, col="royalblue3", font=2)
	#DNAPK
	points(normdata_sorted["P78527", "rownum"], normdata_sorted["P78527","mock_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P78527", "rownum"], normdata_sorted["P78527","mock_avgzscore"], labels="PRKDC", cex=0.75, pos=2, col="royalblue3", font=2)
	#IFI16
	points(normdata_sorted["Q16666", "rownum"], normdata_sorted["Q16666","mock_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q16666", "rownum"], normdata_sorted["Q16666","mock_avgzscore"], labels="IFI16", cex=0.75, pos=2, col="royalblue3", font=2)
	#ATRX
	points(normdata_sorted["P46100", "rownum"], normdata_sorted["P46100","mock_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P46100", "rownum"], normdata_sorted["P46100","mock_avgzscore"], labels="ATRX", cex=0.75, pos=3, col="royalblue3", font=2)
	#USP7
	points(normdata_sorted["Q93009", "rownum"], normdata_sorted["Q93009","mock_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q93009", "rownum"], normdata_sorted["Q93009","mock_avgzscore"], labels="USP7", cex=0.75, pos=4, col="royalblue3", font=2)
	#SLFN5
	points(normdata_sorted["Q08AF3", "rownum"], normdata_sorted["Q08AF3","mock_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q08AF3", "rownum"], normdata_sorted["Q08AF3","mock_avgzscore"], labels="SLFN5", cex=0.75, pos=2, col="royalblue3", font=2)
	#SAMHD1
	points(normdata_sorted["Q9Y3Z3", "rownum"], normdata_sorted["Q9Y3Z3","mock_avgzscore"], pch=20, cex=4, col="red")
	#text(normdata_sorted["Q9Y3Z3", "rownum"], normdata_sorted["Q9Y3Z3","mock_avgzscore"], labels="SAMHD1", cex=0.75, pos=2, col="red", font=2)
	
	#axis(1,col="black",cex.axis=0.75)
	mtext("proteins (ordered by zscore)",1,line=0.5)
	axis(2,c(seq(-8,8,by=1)),col="black",cex.axis=0.75)
	mtext("mock zscore",2,line=2.5)
	dev.off()

	pdf("HSV_SnakePlot_wICP0substrates_20171017.pdf")
	mockzscoredata<-normdata[,"hsv_avgzscore"]
	normdata_sorted<-normdata[order(normdata$hsv_avgzscore, decreasing=TRUE), ]
	normdata_sorted$rownum<-1:nrow(normdata_sorted)
	minzscore<-floor(min(normdata_sorted$hsv_avgzscore, na.rm=TRUE))
	maxzscore<-ceiling(max(normdata_sorted$hsv_avgzscore, na.rm=TRUE))
	#if ( abs(minzscore)>maxzscore ){
	#	maxzscore<-abs(minzscore)
	#}else {
	#	minzscore <- -maxzscore
	#}
	par(bty='l',lwd=3)
	plot(normdata_sorted$hsv_avgzscore,
	     pch=20,
	     cex=1.5,
	     col="grey70",
	     ylim=c(minzscore,maxzscore),
	     xlab="", ylab="", xaxt='n', yaxt='n')
	
	#PML
	points(normdata_sorted["P29590","rownum"], normdata_sorted["P29590","hsv_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P29590","rownum"], normdata_sorted["P29590","hsv_avgzscore"], labels="PML", cex=0.75, pos=4, col="royalblue3", font=2)
	#DNAPK
	points(normdata_sorted["P78527", "rownum"], normdata_sorted["P78527","hsv_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P78527", "rownum"], normdata_sorted["P78527","hsv_avgzscore"], labels="PRKDC", cex=0.75, pos=2, col="royalblue3", font=2)
	#IFI16
	points(normdata_sorted["Q16666", "rownum"], normdata_sorted["Q16666","hsv_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q16666", "rownum"], normdata_sorted["Q16666","hsv_avgzscore"], labels="IFI16", cex=0.75, pos=2, col="royalblue3", font=2)
	#ATRX
	points(normdata_sorted["P46100", "rownum"], normdata_sorted["P46100","hsv_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P46100", "rownum"], normdata_sorted["P46100","hsv_avgzscore"], labels="ATRX", cex=0.75, pos=3, col="royalblue3", font=2)
	#USP7
	points(normdata_sorted["Q93009", "rownum"], normdata_sorted["Q93009","hsv_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q93009", "rownum"], normdata_sorted["Q93009","hsv_avgzscore"], labels="USP7", cex=0.75, pos=4, col="royalblue3", font=2)
	#SLFN5
	points(normdata_sorted["Q08AF3", "rownum"], normdata_sorted["Q08AF3","hsv_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q08AF3", "rownum"], normdata_sorted["Q08AF3","hsv_avgzscore"], labels="SLFN5", cex=0.75, pos=2, col="royalblue3", font=2)
	#SAMHD1
	points(normdata_sorted["Q9Y3Z3", "rownum"], normdata_sorted["Q9Y3Z3","hsv_avgzscore"], pch=20, cex=4, col="red")
	#text(normdata_sorted["Q9Y3Z3", "rownum"], normdata_sorted["Q9Y3Z3","hsv_avgzscore"], labels="SAMHD1", cex=0.75, pos=2, col="red", font=2)
	
	#axis(1,col="black",cex.axis=0.75)
	mtext("proteins (ordered by fold change)",1,line=0.5)
	axis(2,c(seq(-8,8,by=1)),col="black",cex.axis=0.75)
	mtext("mock zscore",2,line=2.5)
	dev.off()
	
	pdf("dICP0_SnakePlot_wICP0substrates_20171017.pdf")
	mockzscoredata<-normdata[,"dICP0_avgzscore"]
	normdata_sorted<-normdata[order(normdata$dICP0_avgzscore, decreasing=TRUE), ]
	normdata_sorted$rownum<-1:nrow(normdata_sorted)
	minzscore<-floor(min(normdata_sorted$dICP0_avgzscore, na.rm=TRUE))
	maxzscore<-ceiling(max(normdata_sorted$dICP0_avgzscore, na.rm=TRUE))
	#if ( abs(minzscore)>maxzscore ){
	#	maxzscore<-abs(minzscore)
	#}else {
	#	minzscore <- -maxzscore
	#}
	par(bty='l',lwd=3)
	plot(normdata_sorted$dICP0_avgzscore,
	     pch=20,
	     cex=1.5,
	     col="grey70",
	     ylim=c(minzscore,maxzscore),
	     xlab="", ylab="", xaxt='n', yaxt='n')
	
	#PML
	points(normdata_sorted["P29590","rownum"], normdata_sorted["P29590","dICP0_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P29590","rownum"], normdata_sorted["P29590","dICP0_avgzscore"], labels="PML", cex=0.75, pos=4, col="royalblue3", font=2)
	#DNAPK
	points(normdata_sorted["P78527", "rownum"], normdata_sorted["P78527","dICP0_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P78527", "rownum"], normdata_sorted["P78527","dICP0_avgzscore"], labels="PRKDC", cex=0.75, pos=2, col="royalblue3", font=2)
	#IFI16
	points(normdata_sorted["Q16666", "rownum"], normdata_sorted["Q16666","dICP0_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q16666", "rownum"], normdata_sorted["Q16666","dICP0_avgzscore"], labels="IFI16", cex=0.75, pos=2, col="royalblue3", font=2)
	#ATRX
	points(normdata_sorted["P46100", "rownum"], normdata_sorted["P46100","dICP0_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["P46100", "rownum"], normdata_sorted["P46100","dICP0_avgzscore"], labels="ATRX", cex=0.75, pos=3, col="royalblue3", font=2)
	#USP7
	points(normdata_sorted["Q93009", "rownum"], normdata_sorted["Q93009","dICP0_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q93009", "rownum"], normdata_sorted["Q93009","dICP0_avgzscore"], labels="USP7", cex=0.75, pos=4, col="royalblue3", font=2)
	#SLFN5
	points(normdata_sorted["Q08AF3", "rownum"], normdata_sorted["Q08AF3","dICP0_avgzscore"], pch=20, cex=4, col="royalblue3")
	#text(normdata_sorted["Q08AF3", "rownum"], normdata_sorted["Q08AF3","dICP0_avgzscore"], labels="SLFN5", cex=0.75, pos=2, col="royalblue3", font=2)
	#SAMHD1
	points(normdata_sorted["Q9Y3Z3", "rownum"], normdata_sorted["Q9Y3Z3","dICP0_avgzscore"], pch=20, cex=4, col="red")
	#text(normdata_sorted["Q9Y3Z3", "rownum"], normdata_sorted["Q9Y3Z3","dICP0_avgzscore"], labels="SAMHD1", cex=0.75, pos=2, col="red", font=2)
	
	#axis(1,col="black",cex.axis=0.75)
	mtext("proteins (ordered by fold change)",1,line=0.5)
	axis(2,c(seq(-8,8,by=1)),col="black",cex.axis=0.75)
	mtext("mock zscore",2,line=2.5)
	dev.off()
	
	return()

}



##### MAIN #####


##directory where analysis files are written
analysisdir<-"analysisoutput_AUG172019"


##read raw data
#formatted with experiments as columns and samples as rows
rawdata <- read.table("./datafiles/iPONDdata_NoImputation/dICP0_iBAQvals_Raw_20200817.txt", sep="\t", row.names=1, header=TRUE)
#protein id to gene name table
geneiddata <- read.table("./datafiles/IDfiles/dICP0_ProteinToGeneName.txt", sep="\t", row.names=1, header=TRUE)
#merge rawdata and geneiddata
rawdata <- merge(rawdata, geneiddata, by="row.names", all.x=TRUE)
row.names(rawdata)=rawdata$Row.names
rawdata$Row.names <- NULL
#protein id to protein name table
proteinnamedata <- read.table("./datafiles/IDfiles/dICP0_ProteinToProteinName.txt", sep="\t", row.names=1, header=TRUE)
#merge rawdata and proteinnamedata
rawdata <- merge(rawdata, proteinnamedata, by="row.names", all.x=TRUE)
row.names(rawdata)=rawdata$Row.names
rawdata$Row.names <- NULL
head(rawdata)
rawdata <- rawdata[c("geneID", "protNAME", "mock_1", "mock_2", "mock_3", "HSV_1", "HSV_2", "HSV_3", "dICP0_1", "dICP0_2", "dICP0_3")]
head(rawdata)
#
##end read Log2 transformed, imputed, normalized quantification data
#
##change all "0" intensities to "NA" to avoid influencing the normalization
rawdata[, "mock_1"][rawdata[, "mock_1"] == 0] <- NA
rawdata[, "mock_2"][rawdata[, "mock_2"] == 0] <- NA
rawdata[, "mock_3"][rawdata[, "mock_3"] == 0] <- NA
rawdata[, "HSV_1"][rawdata[, "HSV_1"] == 0] <- NA
rawdata[, "HSV_2"][rawdata[, "HSV_2"] == 0] <- NA
rawdata[, "HSV_3"][rawdata[, "HSV_3"] == 0] <- NA
rawdata[, "dICP0_1"][rawdata[, "dICP0_1"] == 0] <- NA
rawdata[, "dICP0_2"][rawdata[, "dICP0_2"] == 0] <- NA
rawdata[, "dICP0_3"][rawdata[, "dICP0_3"] == 0] <- NA
head(rawdata)
#
##normalize...
##make boxplot of non-normalized raw values
quantcols <- c('mock_1', 'mock_2', 'mock_3', 'HSV_1', 'HSV_2', 'HSV_3', 'dICP0_1', 'dICP0_2', 'dICP0_3')
rawdataMat <- as.matrix(rawdata[,quantcols])
pdf(file.path(".", analysisdir, "iPOND_RawData_boxplot.pdf"))
boxplot(rawdataMat,las=2,cex.axis=0.5)
dev.off()
#
##transform to log2
for ( i in 1:length(quantcols) ){
  currprotcol <- quantcols[i]
  newprotcol <- paste(currprotcol, '_Log2', sep="")
  rawdata[,newprotcol] <- apply(rawdata, MARGIN=1, FUN=function(x) (log2(as.numeric(x[[currprotcol]]))))
}
##convert any -Inf to NA
rawdata[rawdata=="-Inf"] <- NA
#dim(rawdata)
##print
#head(rawdata)
#
##make boxplot of non-normalized log2 values
log2quantcols <- c('mock_1_Log2', 'mock_2_Log2', 'mock_3_Log2', 'HSV_1_Log2', 'HSV_2_Log2', 'HSV_3_Log2', 'dICP0_1_Log2', 'dICP0_2_Log2', 'dICP0_3_Log2')
rawdataLog2Mat <- as.matrix(rawdata[,log2quantcols])
pdf(file.path(".", analysisdir, "iPOND_Log2Data_boxplot.pdf"))
boxplot(rawdataLog2Mat,las=2,cex.axis=0.5)
dev.off()
#
##normalize NonNormLog2 values by median column value
for ( i in 1:length(log2quantcols) ){
  currlog2protcol <- log2quantcols[i]
  newlog2protcol <- paste(currlog2protcol, '_MedNorm', sep="")
  currprotmed <- median(rawdata[,currlog2protcol],na.rm=TRUE)
  rawdata[,newlog2protcol]<-apply(rawdata, MARGIN=1, FUN=function(x) (as.numeric(x[[currlog2protcol]])-currprotmed))
}
#dim(rawdata)
#
##print
#head(rawdata)
##make boxplot of normalized log2 values
log2normquantcols <- c('mock_1_Log2_MedNorm', 'mock_2_Log2_MedNorm', 'mock_3_Log2_MedNorm', 'HSV_1_Log2_MedNorm', 'HSV_2_Log2_MedNorm', 'HSV_3_Log2_MedNorm', 'dICP0_1_Log2_MedNorm', 'dICP0_2_Log2_MedNorm', 'dICP0_3_Log2_MedNorm')
rawdataLog2NormMat <- as.matrix(rawdata[,log2normquantcols])
pdf(file.path(".", analysisdir, "iPOND_Log2NormData_boxplot.pdf"))
boxplot(rawdataLog2NormMat,las=2,cex.axis=0.5)
dev.off()
#
#head(rawdata)
#dim(rawdata)



##read Log2 transformed, imputed, normalized quantification data
#
#formatted with experiments as columns and samples as rows
imputeddata <- read.table("./datafiles/iPONDdata_Imputation/dICP0_iBAQvals_Raw_WithImputations_Log2NormalizedMedian_20160323.txt", sep="\t", row.names=1, header=TRUE)
#protein id to gene name table
geneiddata <- read.table("./datafiles/IDfiles/dICP0_ProteinToGeneName.txt", sep="\t", row.names=1, header=TRUE)
#merge imputeddata and geneiddata
imputeddata <- merge(imputeddata, geneiddata, by="row.names", all.x=TRUE)
row.names(imputeddata)=imputeddata$Row.names
imputeddata$Row.names <- NULL
#protein id to protein name table
proteinnamedata <- read.table("./datafiles/IDfiles/dICP0_ProteinToProteinName.txt", sep="\t", row.names=1, header=TRUE)
#merge imputeddata and proteinnamedata
imputeddata <- merge(imputeddata, proteinnamedata, by="row.names", all.x=TRUE)
row.names(imputeddata)=imputeddata$Row.names
imputeddata$Row.names <- NULL
#
##end read Log2 transformed, imputed, normalized quantification data


##remove background proteins from quantification dataset
#
#read background proteins data
#mock background
backgroundproteins_mock <- read.table("/Users/dybasj/JoeWork/iPONDbackgroundproteins/iPONDMinusBiotinBackground_Mock.txt", sep="\t", row.names=1, header=FALSE)
#HSV background
backgroundproteins_hsv <- read.table("/Users/dybasj/JoeWork/iPONDbackgroundproteins/iPONDMinusBiotinBackground_HSV1.txt", sep="\t", row.names=1, header=FALSE)
#combine mock and HSV
backgroundproteins <- unique(rbind(backgroundproteins_mock, backgroundproteins_hsv))
backgroundproteinids <- row.names(backgroundproteins)
#remove proteins from raw dataset
rawdata<-rawdata[!row.names(rawdata)%in%backgroundproteinids,]
#remove proteins from dataset
imputeddata<-imputeddata[!row.names(imputeddata)%in%backgroundproteinids,]
#
##end remove background proteins from quantification dataset


##perform pca plot on imputed data
#
print("perform PCA")
PCAdata<-imputeddata[,c("mock_1", "mock_2", "mock_3", "HSV_1", "HSV_2", "HSV_3", "dICP0_1", "dICP0_2", "dICP0_3")]
dim(PCAdata)
PCAdataFiltered<-na.omit(PCAdata)
#head(PCAdata)
#print(PCAdataFiltered)
dim(PCAdataFiltered)
#
#to subset the data
#proteindata<-proteindata[!(proteindata$proteinID=="mock_1" | proteindata$proteinID=="mock_2" | proteindata$proteinID=="mock_3"),]
#print(proteindata)
#dim(proteindata)
#
##separate protein ids and quantification data
##########PCAdata_proteins <- PCAdataFiltered[,1]
PCAdata_proteins <- row.names(PCAdataFiltered)
PCAdata_data <- PCAdataFiltered[,1:dim(PCAdataFiltered)[2]]
#
##perform pca using prcomp
PCAdata_pca <- prcomp(PCAdata_data, center=TRUE, scale=TRUE)
print(PCAdata_pca)
summary(PCAdata_pca)
#
test<-PCAdata_pca[["rotation"]]
#
print(test)
#xmin <- 0.39 
#xmax <- 0.42 
#ymin <- -0.5 
#ymax <- 0.7
pdf("MockWTdICP0_PCAplot.pdf")
par(bty="l")
plot(test[,"PC1"],test[,"PC2"],
     #pch=c(16,18,16,18,16,18),
	 cex=3,
	 col="gray50",
	 #col=c("goldenrod1","goldenrod1","royalblue4","royalblue4","dodgerblue","dodgerblue"),
	 #xlab="",ylab="",
	 #xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 #xaxt='n',yaxt='n'
	)
text(test[,"PC1"], test[,"PC2"], labels=row.names(test), cex=1)
#label and format axes
#axis(1, c(seq(xmin,xmax,by=0.0025)), labels=F, col="black",cex.axis=1, tck=-0.01)
#axis(1, c(seq(xmin,xmax,by=0.005)), labels=T, col="black",cex.axis=1)
mtext("PC1",1,line=2.75,cex=1.25)
#axis(2, c(seq(ymin,ymax,by=0.1)), labels=F, col="black",cex.axis=1, tck=-0.01)
#axis(2, c(seq(ymin,ymax,by=0.2)), labels=T, col="black",cex.axis=1)
mtext("PC2",2,line=2.75,cex=1.25)
dev.off()
#
##end perform pca plot


##raw data count replicates with identified protein
##count replicates by using sum of non-NA intensity values
#mock replicate counts
rawdata$mock_NumReps <- apply(rawdata, MARGIN=1, FUN=function(x) sum(!is.na(x[["mock_1"]]),!is.na(x[["mock_2"]]),!is.na(x[["mock_3"]])))
#HSV replicate counts
rawdata$HSV_NumReps <- apply(rawdata, MARGIN=1, FUN=function(x) sum(!is.na(x[["HSV_1"]]),!is.na(x[["HSV_2"]]),!is.na(x[["HSV_3"]])))
#dICP0 replicate counts
rawdata$dICP0_NumReps <- apply(rawdata, MARGIN=1, FUN=function(x) sum(!is.na(x[["dICP0_1"]]),!is.na(x[["dICP0_2"]]),!is.na(x[["dICP0_3"]])))
#
#get list of all proteins identified by iPOND
rawdata_allid <- subset(rawdata, (rawdata$mock_NumReps>=2 | rawdata$HSV_NumReps>=2 | rawdata$dICP0_NumReps>=2))
dim(rawdata_allid)
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_MockWTdICP0_AllIdentifiedProteins_20200817.txt")
write.table(rawdata_allid, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
head(rawdata)
##end raw data count replicates with identified protein



##raw data average replicate abundance values
#mock
rawdata$mock_Log2MedNorm_Avg<-apply(rawdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["mock_1_Log2_MedNorm"]],x[["mock_2_Log2_MedNorm"]],x[["mock_3_Log2_MedNorm"]])), na.rm=TRUE))
#HSV
rawdata$HSV_Log2MedNorm_Avg<-apply(rawdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["HSV_1_Log2_MedNorm"]],x[["HSV_2_Log2_MedNorm"]],x[["HSV_3_Log2_MedNorm"]])), na.rm=TRUE))
#dICP0
rawdata$dICP0_Log2MedNorm_Avg<-apply(rawdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["dICP0_1_Log2_MedNorm"]],x[["dICP0_2_Log2_MedNorm"]],x[["dICP0_3_Log2_MedNorm"]])), na.rm=TRUE))


##zscores for raw data mock, hsv, dICP0
#
print("zscores")
#mock_1 mean, sd
raw_mock_1_mean<-mean(as.numeric(rawdata$mock_1_Log2_MedNorm), na.rm=TRUE)
raw_mock_1_stdev<-sd(as.numeric(rawdata$mock_1_Log2_MedNorm), na.rm=TRUE)
#raw_mock_1 zscore
rawdata$mock_1_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["mock_1_Log2_MedNorm"]]), raw_mock_1_mean, raw_mock_1_stdev))
#
#mock_2 mean, sd
raw_mock_2_mean<-mean(as.numeric(rawdata$mock_2_Log2_MedNorm), na.rm=TRUE)
raw_mock_2_stdev<-sd(as.numeric(rawdata$mock_2_Log2_MedNorm), na.rm=TRUE)
#mock_2 zscore
rawdata$mock_2_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["mock_2_Log2_MedNorm"]]), raw_mock_2_mean, raw_mock_2_stdev))
#
#mock_3 mean, sd
raw_mock_3_mean<-mean(as.numeric(rawdata$mock_3_Log2_MedNorm), na.rm=TRUE)
raw_mock_3_stdev<-sd(as.numeric(rawdata$mock_3_Log2_MedNorm), na.rm=TRUE)
#mock_3 zscore
rawdata$mock_3_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["mock_3_Log2_MedNorm"]]), raw_mock_3_mean, raw_mock_3_stdev))
#
#HSV_1 mean, sd
raw_HSV_1_mean<-mean(as.numeric(rawdata$HSV_1_Log2_MedNorm), na.rm=TRUE)
raw_HSV_1_stdev<-sd(as.numeric(rawdata$HSV_1_Log2_MedNorm), na.rm=TRUE)
#HSV_1 zscore
rawdata$HSV_1_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["HSV_1_Log2_MedNorm"]]), raw_HSV_1_mean, raw_HSV_1_stdev))
#
#HSV_2 mean, sd
raw_HSV_2_mean<-mean(as.numeric(rawdata$HSV_2_Log2_MedNorm), na.rm=TRUE)
raw_HSV_2_stdev<-sd(as.numeric(rawdata$HSV_2_Log2_MedNorm), na.rm=TRUE)
#HSV_2 zscore
rawdata$HSV_2_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["HSV_2_Log2_MedNorm"]]), raw_HSV_2_mean, raw_HSV_2_stdev))
#
#HSV_3 mean, sd
raw_HSV_3_mean<-mean(as.numeric(rawdata$HSV_3_Log2_MedNorm), na.rm=TRUE)
raw_HSV_3_stdev<-sd(as.numeric(rawdata$HSV_3_Log2_MedNorm), na.rm=TRUE)
#HSV_3 zscore
rawdata$HSV_3_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["HSV_3_Log2_MedNorm"]]), raw_HSV_3_mean, raw_HSV_3_stdev))
#
#dICP0_1 mean, sd
raw_dICP0_1_mean<-mean(as.numeric(rawdata$dICP0_1_Log2_MedNorm), na.rm=TRUE)
raw_dICP0_1_stdev<-sd(as.numeric(rawdata$dICP0_1_Log2_MedNorm), na.rm=TRUE)
#dICP0_1 zscore
rawdata$dICP0_1_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["dICP0_1_Log2_MedNorm"]]), raw_dICP0_1_mean, raw_dICP0_1_stdev))
#
#dICP0_2 mean, sd
raw_dICP0_2_mean<-mean(as.numeric(rawdata$dICP0_2_Log2_MedNorm), na.rm=TRUE)
raw_dICP0_2_stdev<-sd(as.numeric(rawdata$dICP0_2_Log2_MedNorm), na.rm=TRUE)
#dICP0_2 zscore
rawdata$dICP0_2_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["dICP0_2_Log2_MedNorm"]]), raw_dICP0_2_mean, raw_dICP0_2_stdev))
#
#dICP0_3 mean, sd
raw_dICP0_3_mean<-mean(as.numeric(rawdata$dICP0_3_Log2_MedNorm), na.rm=TRUE)
raw_dICP0_3_stdev<-sd(as.numeric(rawdata$dICP0_3_Log2_MedNorm), na.rm=TRUE)
#dICP0_3 zscore
rawdata$dICP0_3_zscore<-apply(rawdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["dICP0_3_Log2_MedNorm"]]), raw_dICP0_3_mean, raw_dICP0_3_stdev))
#
######calculate the average z-scores for the replicates
#####rawdata$mock_avgzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["mock_NumReps"]]>=2, mean(as.numeric(c(x[["mock_1_zscore"]],x[["mock_2_zscore"]],x[["mock_3_zscore"]])), na.rm=TRUE), NA))
#####rawdata$hsv_avgzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["HSV_NumReps"]]>=2, mean(as.numeric(c(x[["HSV_1_zscore"]],x[["HSV_2_zscore"]],x[["HSV_3_zscore"]])), na.rm=TRUE), NA))
#####rawdata$dICP0_avgzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["dICP0_NumReps"]]>=2, mean(as.numeric(c(x[["dICP0_1_zscore"]],x[["dICP0_2_zscore"]],x[["dICP0_3_zscore"]])), na.rm=TRUE), NA))
#calculate the average z-scores for the replicates
rawdata$mock_avgzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["mock_NumReps"]]>=1, mean(as.numeric(c(x[["mock_1_zscore"]],x[["mock_2_zscore"]],x[["mock_3_zscore"]])), na.rm=TRUE), NA))
rawdata$hsv_avgzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["HSV_NumReps"]]>=1, mean(as.numeric(c(x[["HSV_1_zscore"]],x[["HSV_2_zscore"]],x[["HSV_3_zscore"]])), na.rm=TRUE), NA))
rawdata$dICP0_avgzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["dICP0_NumReps"]]>=1, mean(as.numeric(c(x[["dICP0_1_zscore"]],x[["dICP0_2_zscore"]],x[["dICP0_3_zscore"]])), na.rm=TRUE), NA))
#
#impute -4 zscores for each NA value
rawdata$mock_modifiedzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["mock_NumReps"]]==0, as.numeric(-4.00), ifelse(x[["mock_NumReps"]]>=2, as.numeric(x[["mock_avgzscore"]]), NA)))
rawdata$hsv_modifiedzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["HSV_NumReps"]]==0, as.numeric(-4.00), ifelse(x[["HSV_NumReps"]]>=2, as.numeric(x[["hsv_avgzscore"]]), NA)))
rawdata$dICP0_modifiedzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["dICP0_NumReps"]]==0, as.numeric(-4.00), ifelse(x[["dICP0_NumReps"]]>=2, as.numeric(x[["dICP0_avgzscore"]]), NA)))
#####rawdata$mock_modifiedzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["mock_NumReps"]]==0, as.numeric(-4.00), x[["mock_avgzscore"]]))
#####rawdata$hsv_modifiedzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["HSV_NumReps"]]==0, as.numeric(-4.00), x[["hsv_avgzscore"]]))
#####rawdata$dICP0_modifiedzscore<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse(x[["dICP0_NumReps"]]==0, as.numeric(-4.00), x[["dICP0_avgzscore"]]))
#
#####rawdata_TMP <- data.frame(rawdata)
#####rawdata_TMP[, "mock_avgzscore"][is.na(rawdata_TMP[, "mock_avgzscore"])] <- -4
#####rawdata_TMP[, "hsv_avgzscore"][is.na(rawdata_TMP[, "hsv_avgzscore"])] <- -4
#####rawdata_TMP[, "dICP0_avgzscore"][is.na(rawdata_TMP[, "dICP0_avgzscore"])] <- -4
######raw: print all data
#####outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_AllDataWithImputedZscores_20200817.txt")
#####write.table(rawdata_TMP, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#

##end zscores for raw data mock, hsv, dICP0





#pdf(file.path(".", analysisdir, "HSVwtHSVdICP0_RawData_zscores_new_labels.pdf"))
#par(bty='l',lwd=3)
#xmin <- floor(min(as.numeric(rawdata$hsv_modifiedzscore), na.rm=TRUE))
#xmax <- ceiling(max(as.numeric(rawdata$hsv_modifiedzscore), na.rm=TRUE))
#ymin <- floor(min(as.numeric(rawdata$dICP0_modifiedzscore), na.rm=TRUE))
##ymax <- ceiling(max(as.numeric(rawdata$dICP0_modifiedzscore), na.rm=TRUE))
#xmin <- -2
#xmax <- 0
#ymin <- -4.5
#ymax <- -3.5
#plot(rawdata$hsv_modifiedzscore, rawdata$dICP0_modifiedzscore,    
#     pch=20,
#     cex=1.2,
#     col="grey70",
#     xlim=c(xmin,xmax), ylim=c(ymin,ymax),
#     xlab="",ylab="",xaxt='n',yaxt='n'
#)
#text(rawdata$hsv_modifiedzscore,rawdata$dICP0_modifiedzscore, labels=rawdata$geneID, cex=0.5, pos=4, col="red", font=2)
##
###add known substrates
##for ( i in 1:length(ICP0subs) ){
##  currsub<-ICP0subs[i]
##  print(currsub)
##  currgene<-rawdata[currsub,"geneID"]
###  print(currgene)
#  points(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], pch=20, cex=1.25, col="red")
#  #text(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], labels=rawdata[currsub,"geneID"], cex=0.5, pos=2, col="red", font=2)
#}
##
##add zinc finger proteins
#for ( i in 1:length(ZincFingerProts) ){
#  currsub<-ZincFingerProts[i]
#  points(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], pch=20, cex=1.25, col="blue")
#  #text(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], labels=rawdata[currsub,"geneID"], cex=0.5, pos=1, col="blue", font=2)
#}
#
#axis(1,c(seq(xmin,xmax,by=1)),col="black",cex.axis=1,lwd=3)
#mtext("HSV wt z-score",1,line=2.5)
#axis(2,c(seq(ymin,ymax,by=0.5)),col="black",cex.axis=1,lwd=3)
#mtext("HSV dICP0 z-score",2,line=2.5)
#dev.off()













##zscores for imputed data mock, hsv, dICP0
#
print("zscores")
#mock_1 mean, sd
imputed_mock_1_mean<-mean(as.numeric(imputeddata$mock_1), na.rm=TRUE)
imputed_mock_1_stdev<-sd(as.numeric(imputeddata$mock_1), na.rm=TRUE)
#mock_1 zscore
imputeddata$mock_1_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["mock_1"]]), imputed_mock_1_mean, imputed_mock_1_stdev))
#
#mock_2 mean, sd
imputed_mock_2_mean<-mean(as.numeric(imputeddata$mock_2), na.rm=TRUE)
imputed_mock_2_stdev<-sd(as.numeric(imputeddata$mock_2), na.rm=TRUE)
#mock_2 zscore
imputeddata$mock_2_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["mock_2"]]), imputed_mock_2_mean, imputed_mock_2_stdev))
#
#mock_3 mean, sd
imputed_mock_3_mean<-mean(as.numeric(imputeddata$mock_3), na.rm=TRUE)
imputed_mock_3_stdev<-sd(as.numeric(imputeddata$mock_3), na.rm=TRUE)
#mock_3 zscore
imputeddata$mock_3_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["mock_3"]]), imputed_mock_3_mean, imputed_mock_3_stdev))
#
#HSV_1 mean, sd
imputed_HSV_1_mean<-mean(as.numeric(imputeddata$HSV_1), na.rm=TRUE)
imputed_HSV_1_stdev<-sd(as.numeric(imputeddata$HSV_1), na.rm=TRUE)
#HSV_1 zscore
imputeddata$HSV_1_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["HSV_1"]]), imputed_HSV_1_mean, imputed_HSV_1_stdev))
#
#HSV_2 mean, sd
imputed_HSV_2_mean<-mean(as.numeric(imputeddata$HSV_2), na.rm=TRUE)
imputed_HSV_2_stdev<-sd(as.numeric(imputeddata$HSV_2), na.rm=TRUE)
#HSV_2 zscore
imputeddata$HSV_2_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["HSV_2"]]), imputed_HSV_2_mean, imputed_HSV_2_stdev))
#
#HSV_3 mean, sd
imputed_HSV_3_mean<-mean(as.numeric(imputeddata$HSV_3), na.rm=TRUE)
imputed_HSV_3_stdev<-sd(as.numeric(imputeddata$HSV_3), na.rm=TRUE)
#HSV_3 zscore
imputeddata$HSV_3_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["HSV_3"]]), imputed_HSV_3_mean, imputed_HSV_3_stdev))
#
#dICP0_1 mean, sd
imputed_dICP0_1_mean<-mean(as.numeric(imputeddata$dICP0_1), na.rm=TRUE)
imputed_dICP0_1_stdev<-sd(as.numeric(imputeddata$dICP0_1), na.rm=TRUE)
#dICP0_1 zscore
imputeddata$dICP0_1_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["dICP0_1"]]), imputed_dICP0_1_mean, imputed_dICP0_1_stdev))
#
#dICP0_2 mean, sd
imputed_dICP0_2_mean<-mean(as.numeric(imputeddata$dICP0_2), na.rm=TRUE)
imputed_dICP0_2_stdev<-sd(as.numeric(imputeddata$dICP0_2), na.rm=TRUE)
#dICP0_2 zscore
imputeddata$dICP0_2_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["dICP0_2"]]), imputed_dICP0_2_mean, imputed_dICP0_2_stdev))
#
#dICP0_3 mean, sd
imputed_dICP0_3_mean<-mean(as.numeric(imputeddata$dICP0_3), na.rm=TRUE)
imputed_dICP0_3_stdev<-sd(as.numeric(imputeddata$dICP0_3), na.rm=TRUE)
#dICP0_3 zscore
imputeddata$dICP0_3_zscore<-apply(imputeddata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["dICP0_3"]]), imputed_dICP0_3_mean, imputed_dICP0_3_stdev))
#
#calculate the average z-scores for the replicates
imputeddata$mock_avgzscore<-apply(imputeddata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["mock_1_zscore"]],x[["mock_2_zscore"]],x[["mock_3_zscore"]])), na.rm=TRUE))
imputeddata$hsv_avgzscore<-apply(imputeddata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["HSV_1_zscore"]],x[["HSV_2_zscore"]],x[["HSV_3_zscore"]])), na.rm=TRUE))
imputeddata$dICP0_avgzscore<-apply(imputeddata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["dICP0_1_zscore"]],x[["dICP0_2_zscore"]],x[["dICP0_3_zscore"]])), na.rm=TRUE))
#
###end zscores for imputed data mock, hsv, dICP0


##fold changes for raw mock, hsv, dICP0
#
print("raw data fold changes")
fcresults<-apply(rawdata, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["mock_1_Log2_MedNorm"]],x[["mock_2_Log2_MedNorm"]],x[["mock_3_Log2_MedNorm"]])),as.numeric(c(x[["HSV_1_Log2_MedNorm"]],x[["HSV_2_Log2_MedNorm"]],x[["HSV_3_Log2_MedNorm"]]))))
rawdata$mock_hsv_1_fc<-fcresults[1,]
rawdata$mock_hsv_2_fc<-fcresults[2,]
rawdata$mock_hsv_3_fc<-fcresults[3,]
fcresults<-apply(rawdata, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["mock_1_Log2_MedNorm"]],x[["mock_2_Log2_MedNorm"]],x[["mock_3_Log2_MedNorm"]])),as.numeric(c(x[["dICP0_1_Log2_MedNorm"]],x[["dICP0_2_Log2_MedNorm"]],x[["dICP0_3_Log2_MedNorm"]]))))
rawdata$mock_dICP0_1_fc<-fcresults[1,]
rawdata$mock_dICP0_2_fc<-fcresults[2,]
rawdata$mock_dICP0_3_fc<-fcresults[3,]
fcresults<-apply(rawdata, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["HSV_1_Log2_MedNorm"]],x[["HSV_2_Log2_MedNorm"]],x[["HSV_3_Log2_MedNorm"]])),as.numeric(c(x[["dICP0_1_Log2_MedNorm"]],x[["dICP0_2_Log2_MedNorm"]],x[["dICP0_3_Log2_MedNorm"]]))))
rawdata$hsv_dICP0_1_fc<-fcresults[1,]
rawdata$hsv_dICP0_2_fc<-fcresults[2,]
rawdata$hsv_dICP0_3_fc<-fcresults[3,]
#mean fold change for three replicates of mock, hsv, dICP0
###rawdata$mock_hsv_avgfc<-apply(rawdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["mock_hsv_1_fc"]],x[["mock_hsv_2_fc"]],x[["mock_hsv_3_fc"]])), na.rm=TRUE))
###rawdata$mock_dICP0_avgfc<-apply(rawdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["mock_dICP0_1_fc"]],x[["mock_dICP0_2_fc"]],x[["mock_dICP0_3_fc"]])), na.rm=TRUE))
###rawdata$hsv_dICP0_avgfc<-apply(rawdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["hsv_dICP0_1_fc"]],x[["hsv_dICP0_2_fc"]],x[["hsv_dICP0_3_fc"]])), na.rm=TRUE))
rawdata$mock_hsv_avgfc<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse((x[["mock_NumReps"]]>=2 & x[["HSV_NumReps"]]>=2), foldchange(as.numeric(c(x[["mock_Log2MedNorm_Avg"]])),as.numeric(c(x[["HSV_Log2MedNorm_Avg"]]))), NA))
rawdata$mock_dICP0_avgfc<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse((x[["mock_NumReps"]]>=2 & x[["dICP0_NumReps"]]>=2), foldchange(as.numeric(c(x[["mock_Log2MedNorm_Avg"]])),as.numeric(c(x[["dICP0_Log2MedNorm_Avg"]]))), NA))
rawdata$hsv_dICP0_avgfc<-apply(rawdata, MARGIN=1, FUN=function(x) ifelse((x[["HSV_NumReps"]]>=2 & x[["dICP0_NumReps"]]>=2), foldchange(as.numeric(c(x[["HSV_Log2MedNorm_Avg"]])),as.numeric(c(x[["dICP0_Log2MedNorm_Avg"]]))), NA))
#
##end fold changes for raw mock, hsv, dICP0


##fold changes for imputed mock, hsv, dICP0
#
print("imputated data fold changes")
fcresults<-apply(imputeddata, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["mock_1"]],x[["mock_2"]],x[["mock_3"]])),as.numeric(c(x[["HSV_1"]],x[["HSV_2"]],x[["HSV_3"]]))))
imputeddata$mock_hsv_1_fc<-fcresults[1,]
imputeddata$mock_hsv_2_fc<-fcresults[2,]
imputeddata$mock_hsv_3_fc<-fcresults[3,]
fcresults<-apply(imputeddata, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["mock_1"]],x[["mock_2"]],x[["mock_3"]])),as.numeric(c(x[["dICP0_1"]],x[["dICP0_2"]],x[["dICP0_3"]]))))
imputeddata$mock_dICP0_1_fc<-fcresults[1,]
imputeddata$mock_dICP0_2_fc<-fcresults[2,]
imputeddata$mock_dICP0_3_fc<-fcresults[3,]
fcresults<-apply(imputeddata, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["HSV_1"]],x[["HSV_2"]],x[["HSV_3"]])),as.numeric(c(x[["dICP0_1"]],x[["dICP0_2"]],x[["dICP0_3"]]))))
imputeddata$hsv_dICP0_1_fc<-fcresults[1,]
imputeddata$hsv_dICP0_2_fc<-fcresults[2,]
imputeddata$hsv_dICP0_3_fc<-fcresults[3,]
#mean fold change for three replicates of mock, hsv, dICP0
imputeddata$mock_hsv_avgfc<-apply(imputeddata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["mock_hsv_1_fc"]],x[["mock_hsv_2_fc"]],x[["mock_hsv_3_fc"]])), na.rm=TRUE))
imputeddata$mock_dICP0_avgfc<-apply(imputeddata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["mock_dICP0_1_fc"]],x[["mock_dICP0_2_fc"]],x[["mock_dICP0_3_fc"]])), na.rm=TRUE))
imputeddata$hsv_dICP0_avgfc<-apply(imputeddata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["hsv_dICP0_1_fc"]],x[["hsv_dICP0_2_fc"]],x[["hsv_dICP0_3_fc"]])), na.rm=TRUE))
#
##end fold changes for imputed mock, hsv, dICP0


##significance test comparisons on raw mock, hsv, dICP0
#
print("rawdata sig test")
sigtestresults<-apply(rawdata, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["mock_1_Log2_MedNorm"]],x[["mock_2_Log2_MedNorm"]],x[["mock_3_Log2_MedNorm"]])),as.numeric(c(x[["HSV_1_Log2_MedNorm"]],x[["HSV_2_Log2_MedNorm"]],x[["HSV_3_Log2_MedNorm"]])),'UNPAIRED'))
rawdata$mock_hsv_FtestPval<-sigtestresults[1,]
rawdata$mock_hsv_TtestPval<-sigtestresults[2,]
sigtestresults<-apply(rawdata, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["mock_1_Log2_MedNorm"]],x[["mock_2_Log2_MedNorm"]],x[["mock_3_Log2_MedNorm"]])),as.numeric(c(x[["dICP0_1_Log2_MedNorm"]],x[["dICP0_2_Log2_MedNorm"]],x[["dICP0_3_Log2_MedNorm"]])),'UNPAIRED'))
rawdata$mock_dICP0_FtestPval<-sigtestresults[1,]
rawdata$mock_dICP0_TtestPval<-sigtestresults[2,]
sigtestresults<-apply(rawdata, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["HSV_1_Log2_MedNorm"]],x[["HSV_2_Log2_MedNorm"]],x[["HSV_3_Log2_MedNorm"]])),as.numeric(c(x[["dICP0_1_Log2_MedNorm"]],x[["dICP0_2_Log2_MedNorm"]],x[["dICP0_3_Log2_MedNorm"]])),'UNPAIRED'))
rawdata$hsv_dICP0_FtestPval<-sigtestresults[1,]
rawdata$hsv_dICP0_TtestPval<-sigtestresults[2,]
#neg log10 pvalue for mock, hsv and dICP0 comparisons
rawdata$mock_hsv_neglog10pval<-apply(rawdata, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["mock_hsv_TtestPval"]])))
rawdata$mock_dICP0_neglog10pval<-apply(rawdata, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["mock_dICP0_TtestPval"]])))
rawdata$hsv_dICP0_neglog10pval<-apply(rawdata, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["hsv_dICP0_TtestPval"]])))
#
##end significance test comparisons on raw mock, hsv, dICP0


##significance test comparisons on imputed mock, hsv, dICP0
#
print("imputed sig test")
sigtestresults<-apply(imputeddata, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["mock_1"]],x[["mock_2"]],x[["mock_3"]])),as.numeric(c(x[["HSV_1"]],x[["HSV_2"]],x[["HSV_3"]])),'UNPAIRED'))
imputeddata$mock_hsv_FtestPval<-sigtestresults[1,]
imputeddata$mock_hsv_TtestPval<-sigtestresults[2,]
sigtestresults<-apply(imputeddata, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["mock_1"]],x[["mock_2"]],x[["mock_3"]])),as.numeric(c(x[["dICP0_1"]],x[["dICP0_2"]],x[["dICP0_3"]])),'UNPAIRED'))
imputeddata$mock_dICP0_FtestPval<-sigtestresults[1,]
imputeddata$mock_dICP0_TtestPval<-sigtestresults[2,]
sigtestresults<-apply(imputeddata, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["HSV_1"]],x[["HSV_2"]],x[["HSV_3"]])),as.numeric(c(x[["dICP0_1"]],x[["dICP0_2"]],x[["dICP0_3"]])),'UNPAIRED'))
imputeddata$hsv_dICP0_FtestPval<-sigtestresults[1,]
imputeddata$hsv_dICP0_TtestPval<-sigtestresults[2,]
#neg log10 pvalue for mock, hsv and dICP0 comparisons
imputeddata$mock_hsv_neglog10pval<-apply(imputeddata, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["mock_hsv_TtestPval"]])))
imputeddata$mock_dICP0_neglog10pval<-apply(imputeddata, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["mock_dICP0_TtestPval"]])))
imputeddata$hsv_dICP0_neglog10pval<-apply(imputeddata, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["hsv_dICP0_TtestPval"]])))
#
##end significance test comparisons on imputed mock, hsv, dICP0


##set analysis for raw data
#
#proteins enriched on mock compared to hsv
mock_enriched_hsv_prots_raw_data <- subset(rawdata, (((mock_NumReps>=2 & HSV_NumReps>=2)  & (rawdata$mock_hsv_avgfc <= -1 & rawdata$mock_hsv_TtestPval < 0.05)) | (mock_NumReps>=2 & HSV_NumReps==0)))
mock_enriched_hsv_prots_raw <- row.names(mock_enriched_hsv_prots_raw_data)
mock_enriched_hsv_prots_raw_num <- length(mock_enriched_hsv_prots_raw)
print("raw: mock enriched vs hsv")
print(mock_enriched_hsv_prots_raw_num)
#proteins enriched on hsv compared to mock
hsv_enriched_mock_prots_raw_data <- subset(rawdata, (((mock_NumReps>=2 & HSV_NumReps>=2)  & (rawdata$mock_hsv_avgfc >= 1 & rawdata$mock_hsv_TtestPval < 0.05)) | (mock_NumReps==0 & HSV_NumReps>=2)))
hsv_enriched_mock_prots_raw <- row.names(hsv_enriched_mock_prots_raw_data)
hsv_enriched_mock_prots_raw_num <- length(hsv_enriched_mock_prots_raw)
print("raw: hsv enriched vs mock")
print(hsv_enriched_mock_prots_raw_num)
print(hsv_enriched_mock_prots_raw)
#proteins not enriched for mock-hsv comparison
mock_hsv_nonenriched_prots_raw_data <- subset(rawdata, ((mock_NumReps>=2 & HSV_NumReps>=2)  & ((rawdata$mock_hsv_avgfc < 1 & rawdata$mock_hsv_avgfc > -1) | (rawdata$mock_hsv_TtestPval >= 0.05))))
mock_hsv_nonenriched_prots_raw <- row.names(mock_hsv_nonenriched_prots_raw_data)
mock_hsv_nonenriched_prots_raw_num <- length(mock_hsv_nonenriched_prots_raw)
print("raw: mock vs hsv not-enriched")
print(mock_hsv_nonenriched_prots_raw_num)
#
#proteins enriched on mock compared to dICP0
mock_enriched_dICP0_prots_raw_data <- subset(rawdata, (((mock_NumReps>=2 & dICP0_NumReps>=2)  & (rawdata$mock_dICP0_avgfc <= -1 & rawdata$mock_dICP0_TtestPval < 0.05)) | (mock_NumReps>=2 & dICP0_NumReps==0)))
mock_enriched_dICP0_prots_raw <- row.names(mock_enriched_dICP0_prots_raw_data)
mock_enriched_dICP0_prots_raw_num <- length(mock_enriched_dICP0_prots_raw)
print("raw: mock enriched vs dICP0")
print(mock_enriched_dICP0_prots_raw_num)
#proteins enriched on dICP0 compared to mock
dICP0_enriched_mock_prots_raw_data <- subset(rawdata, (((mock_NumReps>=2 & dICP0_NumReps>=2)  & (rawdata$mock_dICP0_avgfc >= 1 & rawdata$mock_dICP0_TtestPval < 0.05)) | (mock_NumReps==0 & dICP0_NumReps>=2)))
dICP0_enriched_mock_prots_raw <- row.names(dICP0_enriched_mock_prots_raw_data)
dICP0_enriched_mock_prots_raw_num <- length(dICP0_enriched_mock_prots_raw)
print("raw: dICP0 enriched vs mock")
print(dICP0_enriched_mock_prots_raw_num)
print(dICP0_enriched_mock_prots_raw)
#proteins not enriched for mock-hsv comparison
mock_dICP0_nonenriched_prots_raw_data <- subset(rawdata, ((mock_NumReps>=2 & dICP0_NumReps>=2)  & ((rawdata$mock_dICP0_avgfc < 1 & rawdata$mock_dICP0_avgfc > -1) | (rawdata$mock_dICP0_TtestPval >= 0.05))))
mock_dICP0_nonenriched_prots_raw <- row.names(mock_dICP0_nonenriched_prots_raw_data)
mock_dICP0_nonenriched_prots_raw_num <- length(mock_dICP0_nonenriched_prots_raw)
print("raw: mock vs dICP0 not-enriched")
print(mock_dICP0_nonenriched_prots_raw_num)
#
#proteins enriched on hsv compared to dICP0
hsv_enriched_dICP0_prots_raw_data <- subset(rawdata, (((HSV_NumReps>=2 & dICP0_NumReps>=2)  & (rawdata$hsv_dICP0_avgfc <= 0 & rawdata$hsv_dICP0_TtestPval < 0.05)) | (HSV_NumReps>=2 & dICP0_NumReps==0)))
hsv_enriched_dICP0_prots_raw <- row.names(hsv_enriched_dICP0_prots_raw_data)
hsv_enriched_dICP0_prots_raw_num <- length(hsv_enriched_dICP0_prots_raw)
print("raw: hsv enriched vs dICP0")
print(hsv_enriched_dICP0_prots_raw_num)
print(hsv_enriched_dICP0_prots_raw)
#proteins enriched on dICP0 compared to hsv
dICP0_enriched_hsv_prots_raw_data <- subset(rawdata, (((HSV_NumReps>=2 & dICP0_NumReps>=2)  & (rawdata$hsv_dICP0_avgfc >= 0 & rawdata$hsv_dICP0_TtestPval < 0.05)) | (HSV_NumReps==0 & dICP0_NumReps>=2)))
dICP0_enriched_hsv_prots_raw <- row.names(dICP0_enriched_hsv_prots_raw_data)
dICP0_enriched_hsv_prots_raw_num <- length(dICP0_enriched_hsv_prots_raw)
print("raw: dICP0 enriched vs hsv")
print(dICP0_enriched_hsv_prots_raw_num)
print(dICP0_enriched_hsv_prots_raw)
#proteins not enriched for hsv-dICP0 comparison
hsv_dICP0_nonenriched_prots_raw_data <- subset(rawdata, ((HSV_NumReps>=2 & dICP0_NumReps>=2)  & ((rawdata$hsv_dICP0_avgfc < 0 & rawdata$hsv_dICP0_avgfc > 0) | (rawdata$hsv_dICP0_TtestPval >= 0.05))))
hsv_dICP0_nonenriched_prots_raw <- row.names(hsv_dICP0_nonenriched_prots_raw_data)
hsv_dICP0_nonenriched_prots_raw_num <- length(hsv_dICP0_nonenriched_prots_raw)
print("raw: hsv vs dICP0 not-enriched")
print(hsv_dICP0_nonenriched_prots_raw_num)
#
#proteins enriched on HSV WT and/or dICP0 compared to mock
hsv_dICP0_enriched_mock_prots_raw <- intersect(hsv_enriched_mock_prots_raw, dICP0_enriched_mock_prots_raw)
hsv_dICP0_enriched_mock_prots_raw_num <- length(hsv_dICP0_enriched_mock_prots_raw)
print("raw: HSV and dICP0 enriched vs mock")
print(hsv_dICP0_enriched_mock_prots_raw_num)
print(hsv_dICP0_enriched_mock_prots_raw)
hsv_enriched_dICP0_notenriched_mock_prots_raw <- setdiff(hsv_enriched_mock_prots_raw, dICP0_enriched_mock_prots_raw)
hsv_enriched_dICP0_notenriched_mock_prots_raw_num <- length(hsv_enriched_dICP0_notenriched_mock_prots_raw)
hsv_enriched_dICP0_notenriched_mock_prots_raw_data <- rawdata[row.names(rawdata)%in%hsv_enriched_dICP0_notenriched_mock_prots_raw,]
print("raw: HSV enriched vs mock, dICP0 not enriched vs mock")
print(hsv_enriched_dICP0_notenriched_mock_prots_raw_num)
print(hsv_enriched_dICP0_notenriched_mock_prots_raw)
dICP0_enriched_hsv_notenriched_mock_prots_raw <- setdiff(dICP0_enriched_mock_prots_raw, hsv_enriched_mock_prots_raw)
dICP0_enriched_hsv_notenriched_mock_prots_raw_num <- length(dICP0_enriched_hsv_notenriched_mock_prots_raw)
dICP0_enriched_hsv_notenriched_mock_prots_raw_data <- rawdata[row.names(rawdata)%in%dICP0_enriched_hsv_notenriched_mock_prots_raw,]
print("raw: dICP0 enriched vs mock, HSV not enriched vs mock")
print(dICP0_enriched_hsv_notenriched_mock_prots_raw_num)
print(dICP0_enriched_hsv_notenriched_mock_prots_raw)
#
#proteins with highest abundance on HSV WT based on HSV z-score
rawdata_HSV_2reps <- subset(rawdata, rawdata$HSV_NumReps>=2)
dim(rawdata_HSV_2reps)
rawdata_HSV_2reps_sorted <- rawdata_HSV_2reps[order(rawdata_HSV_2reps$hsv_avgzscore, decreasing=TRUE),]
rawdata_HSV_highabundance_data <- rawdata_HSV_2reps_sorted[1:650,]
rawdata_HSV_highabundance_prot <- row.names(rawdata_HSV_highabundance_data)
rawdata_HSV_highabundance_gene <- rawdata_HSV_highabundance_data$geneID
rawdata_HSV_highabundance_prot_num <- length(rawdata_HSV_highabundance_prot)
print(rawdata_HSV_highabundance_prot_num)
print(rawdata_HSV_highabundance_gene)
#proteins with highest abundance on dICP0 based on dICP0 z-score
rawdata_dICP0_2reps <- subset(rawdata, rawdata$dICP0_NumReps>=2)
dim(rawdata_dICP0_2reps)
rawdata_dICP0_2reps_sorted <- rawdata_dICP0_2reps[order(rawdata_dICP0_2reps$dICP0_avgzscore, decreasing=TRUE),]
rawdata_dICP0_highabundance_data <- rawdata_dICP0_2reps_sorted[1:650,]
rawdata_dICP0_highabundance_prot <- row.names(rawdata_dICP0_highabundance_data)
rawdata_dICP0_highabundance_gene <- rawdata_dICP0_highabundance_data$geneID
rawdata_dICP0_highabundance_prot_num <- length(rawdata_dICP0_highabundance_prot)
print(rawdata_dICP0_highabundance_prot_num)
print(rawdata_dICP0_highabundance_gene)
rawdata_highabundance_ICP0sub_inhibit<-setdiff(rawdata_dICP0_highabundance_prot,rawdata_HSV_highabundance_prot)
rawdata_highabundance_ICP0sub_inhibit_data<-rawdata[row.names(rawdata)%in%rawdata_highabundance_ICP0sub_inhibit,]
print(rawdata_highabundance_ICP0sub_inhibit)
dim(rawdata_highabundance_ICP0sub_inhibit_data)
rawdata_highabundance_ICP0sub_recruit<-setdiff(rawdata_HSV_highabundance_prot,rawdata_dICP0_highabundance_prot)
rawdata_highabundance_ICP0sub_recruit_data<-rawdata[row.names(rawdata)%in%rawdata_highabundance_ICP0sub_recruit,]
print(rawdata_highabundance_ICP0sub_recruit)
dim(rawdata_highabundance_ICP0sub_recruit_data)
rawdata_highabundance_HSVdCIP0intersect<-intersect(rawdata_HSV_highabundance_prot,rawdata_dICP0_highabundance_prot)
rawdata_highabundance_HSVdCIP0intersect_num<-length(rawdata_highabundance_HSVdCIP0intersect)
#
#proteins unique to HSV WT compared to dICP0
hsv_unique_dICP0_prots_raw_data <- subset(rawdata, (HSV_NumReps>=2 & dICP0_NumReps==0))
hsv_unique_dICP0_prots_raw <- row.names(hsv_unique_dICP0_prots_raw_data)
hsv_unique_dICP0_prots_raw_num <- length(hsv_unique_dICP0_prots_raw)
print("raw: hsv unique vs dICP0")
print(hsv_unique_dICP0_prots_raw_num)
print(hsv_unique_dICP0_prots_raw)
#
#proteins unique to HSV WT compared to dICP0
dICP0_unique_hsv_prots_raw_data <- subset(rawdata, (HSV_NumReps==0 & dICP0_NumReps>=2))
dICP0_unique_hsv_prots_raw <- row.names(dICP0_unique_hsv_prots_raw_data)
dICP0_unique_hsv_prots_raw_num <- length(dICP0_unique_hsv_prots_raw)
print("raw: dICP0 unique vs hsv")
print(dICP0_unique_hsv_prots_raw_num)
print(dICP0_unique_hsv_prots_raw)
#
#proteins that are high abundant or unique on dICP0 compared to HSV WT
ICP0predSub_decreased <- union(rawdata_highabundance_ICP0sub_inhibit, dICP0_unique_hsv_prots_raw)
print("ICP0predSub_decreased")
print(ICP0predSub_decreased)
#
#proteins that are high abundant or unique on HSV WT compared to HSV dICP0
ICP0predSub_recruited <- union(rawdata_highabundance_ICP0sub_recruit, hsv_unique_dICP0_prots_raw)
print("ICP0predSub_recruited")
print(ICP0predSub_recruited)


testenrichedmissed <- setdiff(dICP0_enriched_hsv_prots_raw, ICP0predSub_decreased)
print(testenrichedmissed)

testenrichedmissedrecruit <- setdiff(hsv_enriched_dICP0_prots_raw, ICP0predSub_recruited)
print(testenrichedmissedrecruit)



#test <- Reduce(intersect, list(dICP0_enriched_hsv_prots_raw, dICP0_enriched_hsv_notenriched_mock_prots_raw, rawdata_highabundance_ICP0sub_inhibit))
#print(test)
#test<-intersect(rawdata_highabundance_ICP0sub_inhibit,dICP0_enriched_hsv_prots_raw)
#print(test)
#testII<-intersect(dICP0_enriched_hsv_prots_raw,dICP0_enriched_hsv_notenriched_mock_prots_raw)
#print(testII)
#testIII<-intersect(rawdata_highabundance_ICP0sub_inhibit,dICP0_enriched_hsv_notenriched_mock_prots_raw)
#print(testIII)
##end set analysis for raw data


##set analysis for imputed data
#
#proteins enriched on mock compared to hsv
mock_enriched_hsv_prots_data <- subset(imputeddata, imputeddata$mock_hsv_avgfc <= -1 & imputeddata$mock_hsv_TtestPval < 0.05)
mock_enriched_hsv_prots <- row.names(mock_enriched_hsv_prots_data)
mock_enriched_hsv_prots_num <- length(mock_enriched_hsv_prots)
print("mock enriched vs hsv")
print(mock_enriched_hsv_prots_num)
#proteins enriched on hsv compared to mock
hsv_enriched_mock_prots_data <- subset(imputeddata, imputeddata$mock_hsv_avgfc >= 1 & imputeddata$mock_hsv_TtestPval < 0.05)
hsv_enriched_mock_prots <- row.names(hsv_enriched_mock_prots_data)
hsv_enriched_mock_prots_num <- length(hsv_enriched_mock_prots)
print("hsv enriched vs mock")
print(hsv_enriched_mock_prots_num)
#proteins not enriched for mock-hsv comparison
mock_hsv_nonenriched_prots_data <- subset(imputeddata, (imputeddata$mock_hsv_TtestPval >= 0.05 | (imputeddata$mock_hsv_avgfc > -1 & imputeddata$mock_hsv_avgfc < 1)))
mock_hsv_nonenriched_prots <- row.names(mock_hsv_nonenriched_prots_data)
mock_hsv_nonenriched_prots_num <- length(mock_hsv_nonenriched_prots)
print("mock vs hsv not-enriched")
print(mock_hsv_nonenriched_prots_num)
#
#proteins enriched on mock compared to dICP0
mock_enriched_dICP0_prots_data <- subset(imputeddata, imputeddata$mock_dICP0_avgfc <= -1 & imputeddata$mock_dICP0_TtestPval < 0.05)
mock_enriched_dICP0_prots <- row.names(mock_enriched_dICP0_prots_data)
mock_enriched_dICP0_prots_num <- length(mock_enriched_dICP0_prots)
print("mock enriched vs dICP0")
print(mock_enriched_dICP0_prots_num)
#proteins enriched on dICP0 compared to mock
dICP0_enriched_mock_prots_data <- subset(imputeddata, imputeddata$mock_dICP0_avgfc >= 1 & imputeddata$mock_dICP0_TtestPval < 0.05)
dICP0_enriched_mock_prots <- row.names(dICP0_enriched_mock_prots_data)
dICP0_enriched_mock_prots_num <- length(dICP0_enriched_mock_prots)
print("dICP0 enriched vs mock")
print(dICP0_enriched_mock_prots_num)
#proteins not enriched for mock-dICP0 comparison
mock_dICP0_nonenriched_prots_data <- subset(imputeddata, (imputeddata$mock_dICP0_TtestPval>=0.05 | (imputeddata$mock_dICP0_avgfc > -1 & imputeddata$mock_dICP0_avgfc < 1)))
mock_dICP0_nonenriched_prots <- row.names(mock_dICP0_nonenriched_prots_data)
mock_dICP0_nonenriched_prots_num <- length(mock_dICP0_nonenriched_prots)
print("mock vs dICP0 not-enriched")
print(mock_dICP0_nonenriched_prots_num)
#
#proteins enriched on hsv compared to dICP0
hsv_enriched_dICP0_prots_data <- subset(imputeddata, imputeddata$hsv_dICP0_avgfc < -1 & imputeddata$hsv_dICP0_TtestPval < 0.05)
hsv_enriched_dICP0_prots <- row.names(hsv_enriched_dICP0_prots_data)
hsv_enriched_dICP0_prots_num <- length(hsv_enriched_dICP0_prots)
print("hsv enriched vs dICP0")
print(hsv_enriched_dICP0_prots_num)
print(hsv_enriched_dICP0_prots)
#proteins enriched on dICP0 compared to hsv
dICP0_enriched_hsv_prots_data <- subset(imputeddata, imputeddata$hsv_dICP0_avgfc > 1 & imputeddata$hsv_dICP0_TtestPval < 0.05)
dICP0_enriched_hsv_prots <- row.names(dICP0_enriched_hsv_prots_data)
dICP0_enriched_hsv_prots_num <- length(dICP0_enriched_hsv_prots)
print("dICP0 enriched vs hsv")
print(dICP0_enriched_hsv_prots_num)
print(dICP0_enriched_hsv_prots)
#proteins not enriched for hsv-dICP0 comparison
hsv_dICP0_nonenriched_prots_data <- subset(imputeddata, (imputeddata$hsv_dICP0_TtestPval>=0.05 | (imputeddata$hsv_dICP0_avgfc > -1 & imputeddata$hsv_dICP0_avgfc < 1)))
hsv_dICP0_nonenriched_prots <- row.names(hsv_dICP0_nonenriched_prots_data)
hsv_dICP0_nonenriched_prots_num <- length(hsv_dICP0_nonenriched_prots)
print("hsv vs dICP0 not-enriched")
print(hsv_dICP0_nonenriched_prots_num)
#
##end set analysis


PCAsubClusUnion <-scan("./analysisoutput_AUG2019/PCAclustering/UnionClus_PML-PRKDC-ATRX-IFI16_20200810.txt", character(), quote = "")


##write output
#
#raw: print all data
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_AllData_20200817.txt")
write.table(rawdata, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#raw: print all proteins that are identified on either HSV WT or dICP0 genomes
HSVwtdIPC0identified_data <- subset(rawdata, (rawdata$HSV_NumReps>=2 | rawdata$dICP0_NumReps>=2))
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_HSVwtIDHSVdICP0ID_20200817.txt")
write.table(HSVwtdIPC0identified_data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#raw: print high abundant proteins on HSV dICP0 and not on WT 
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_dICP0highabundance_20200817.txt")
write.table(rawdata_highabundance_ICP0sub_inhibit_data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#raw: print high abundant proteins on HSV WT and not on dICP0 
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_HSVWThighabundance_20200817.txt")
write.table(rawdata_highabundance_ICP0sub_recruit_data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#raw: print proteins enriched on HSV WT compared to mock
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_HSVWTenrichedVsMock-dICP0notenriched_20200817.txt")
write.table(hsv_enriched_dICP0_notenriched_mock_prots_raw_data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#raw: print proteins enriched on dICP0 compared to mock
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_dICP0enrichedVsMock-HSVWTnotenriched_20200817.txt")
write.table(dICP0_enriched_hsv_notenriched_mock_prots_raw_data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#raw: print proteins enriched on HSV WT compared to dICP0
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_HSVWTenrichedVsdICP0_20200817.txt")
write.table(hsv_enriched_dICP0_prots_raw_data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#raw: print proteins enriched on dICP0 compared to HSV WT
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_dICP0enrichedVsHSVWT_20200817.txt")
write.table(dICP0_enriched_hsv_prots_raw_data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#imputed: print all data
#this is the output that contains the data for the supplemental table of SLFN Nat Micro Kim et. al.
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_ImputedDataAnalysis_AllData_20200817.txt")
write.table(imputeddata, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#imputed: print proteins in the union of known ICP0 substrates
imputeddata_ICP0clus <- imputeddata[row.names(imputeddata)%in%PCAsubClusUnion,]
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_ImputedDataAnalysis_ICP0clusUnion_20200817.txt")
write.table(imputeddata_ICP0clus, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
##end write output



##### plotting #####


##venn diagrams
#
#these aren't really venn diagrams. these are the "Weitzman lab" original way of depicting iPOND enriched vs non-enriched
#mock-hsv
print("Mock-HSV Venn")
xnum<-mock_enriched_hsv_prots_num+mock_hsv_nonenriched_prots_num
ynum<-hsv_enriched_mock_prots_num+mock_hsv_nonenriched_prots_num
intnum<-mock_hsv_nonenriched_prots_num
outfile=file.path(".", analysisdir, "Mock-HSV_VennDiagram_20181211.pdf")
plotpairwisevenn(xnum, ynum, intnum, "", "", "gold3", "steelblue4", "gold", "steelblue1", TRUE, outfile)
#mock-dICP0
print("Mock-dICP0 Venn")
xnum<-mock_enriched_dICP0_prots_num+mock_dICP0_nonenriched_prots_num
ynum<-dICP0_enriched_mock_prots_num+mock_dICP0_nonenriched_prots_num
intnum<-mock_dICP0_nonenriched_prots_num
outfile=file.path(".", analysisdir, "Mock-dICP0_VennDiagram_20181211.pdf")
plotpairwisevenn(xnum, ynum, intnum, "", "", "gold3", "darkred", "gold", "lightsalmon", TRUE, outfile)
#hsv-dICP0
print("HSV-dICP0 Venn")
xnum<-hsv_enriched_dICP0_prots_num+hsv_dICP0_nonenriched_prots_num
ynum<-dICP0_enriched_hsv_prots_num+hsv_dICP0_nonenriched_prots_num
intnum<-hsv_dICP0_nonenriched_prots_num
outfile=file.path(".", analysisdir, "HSV-dICP0_VennDiagram_20181211.pdf")
plotpairwisevenn(xnum, ynum, intnum, "", "", "steelblue4", "darkred", "steelblue1", "lightsalmon", TRUE, outfile)
#
##end venn diagrams


##venn diagram comparing HSV WT and dICP0 enriched compared to mock
print("Mock-HSV Venn")
xnum<-hsv_enriched_mock_prots_raw_num
ynum<-dICP0_enriched_mock_prots_raw_num
intnum<-hsv_dICP0_enriched_mock_prots_raw_num
outfile=file.path(".", analysisdir, "HSVWT-dICP0_enrichedmock_VennDiagram_20180817.pdf")
plotpairwisevenn(xnum, ynum, intnum, "", "", "purple", "olivedrab", "purple", "olivedrab", TRUE, outfile)


##venn diagram comparing high abundance proteins on HSV WT and dICP0
print("HSV WT-dICP0 high abundance proteins Venn")
xnum<-rawdata_HSV_highabundance_prot_num
ynum<-rawdata_dICP0_highabundance_prot_num
intnum<-rawdata_highabundance_HSVdCIP0intersect_num
outfile=file.path(".", analysisdir, "HSVWT-dICP0_highabundanceproteins_VennDiagram.pdf")
plotpairwisevenn(xnum, ynum, intnum, "", "", "red", "olivedrab", "red", "olivedrab", TRUE, outfile)


ICP0subs<-c("P29590", "P78527", "Q16666", "P46100", "Q93009", "Q08AF3", "Q14149")

ZincFingerProts<-scan("./iPOND_MockHSVwtHSVdICP0_ZincFingerProteins.txt", character(), quote = "")
ZincFingerProts_Data <- rawdata[row.names(rawdata)%in%ZincFingerProts, ]
dim(ZincFingerProts_Data)
ZincFingerProtsExpessed_Data <- subset(ZincFingerProts_Data, (ZincFingerProts_Data$HSV_NumReps>=2 | ZincFingerProts_Data$dICP0_NumReps>=2))
dim(ZincFingerProtsExpessed_Data)
outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_ZincFinger_20200817.txt")
write.table(ZincFingerProtsExpessed_Data, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)



##########raw data snake plots
#
pdf(file.path(".", analysisdir, "Mock_RawData_SnakePlot.pdf"), height=6, width=8)
rawdata_sorted<-rawdata[order(rawdata$mock_avgzscore, decreasing=TRUE), ]
nazscores <- apply(rawdata_sorted, 1, function(x) is.na(x[["mock_avgzscore"]]))
rawdata_sorted<-rawdata_sorted[ !nazscores, ]
rawdata_sorted$rownum<-1:nrow(rawdata_sorted)
minzscore<-floor(min(rawdata_sorted$mock_avgzscore, na.rm=TRUE))
maxzscore<-ceiling(max(rawdata_sorted$mock_avgzscore, na.rm=TRUE))
#minprotnum<-floor(min(rawdata_sorted$rownum, na.rm=TRUE))
minprotnum<-0
maxprotnum<-ceiling(max(rawdata_sorted$rownum, na.rm=TRUE))
par(bty='l',lwd=3)
plot(rawdata_sorted$mock_avgzscore,
     pch=20,
     cex=1.5,
     col="grey50",
     ylim=c(minzscore,maxzscore),
     xlab="", ylab="", xaxt='n', yaxt='n'
)
#
#color first 650 high abundance proteins
#for ( i in 1:650 ){
#  points(i,rawdata_sorted[i,"mock_avgzscore"], pch=20, cex=1.5, col="skyblue")
#}
#
#add known substrates
for ( i in 1:length(ICP0subs) ){
  currsub<-ICP0subs[i]
  points(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"mock_avgzscore"], pch=20, cex=1.25, col="red")
  #text(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"mock_avgzscore"], labels=rawdata_sorted[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
}
#
#add zinc finger proteins
#for ( i in 1:length(ZincFingerProts) ){
#  currsub<-ZincFingerProts[i]
#  points(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"mock_avgzscore"], pch=20, cex=1.25, col="red")
#  text(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"mock_avgzscore"], labels=rawdata_sorted[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
#}
#
axis(1,c(seq(minprotnum,maxprotnum,by=100)),col="black",cex.axis=1,lwd=3)
mtext("proteins (ordered by zscore)",1,line=2.5)
axis(2,c(seq(minzscore,maxzscore,by=0.5)),col="black",cex.axis=1,lwd=3)
mtext("mock zscore",2,line=2.5)
dev.off()
#
pdf(file.path(".", analysisdir, "HSV_RawData_SnakePlot.pdf"), height=6, width=8)
rawdata_sorted<-rawdata[order(rawdata$hsv_avgzscore, decreasing=TRUE), ]
nazscores <- apply(rawdata_sorted, 1, function(x) is.na(x[["hsv_avgzscore"]]))
rawdata_sorted<-rawdata_sorted[ !nazscores, ]
rawdata_sorted$rownum<-1:nrow(rawdata_sorted)
minzscore<-floor(min(rawdata_sorted$hsv_avgzscore, na.rm=TRUE))
maxzscore<-ceiling(max(rawdata_sorted$hsv_avgzscore, na.rm=TRUE))
#minprotnum<-floor(min(rawdata_sorted$rownum, na.rm=TRUE))
minprotnum<-0
maxprotnum<-ceiling(max(rawdata_sorted$rownum, na.rm=TRUE))
par(bty='l',lwd=3)
plot(rawdata_sorted$hsv_avgzscore,
     pch=20,
     cex=1.5,
     col="grey50",
     ylim=c(minzscore,maxzscore),
     xlab="", ylab="", xaxt='n', yaxt='n'
)
#
#color first 650 high abundance proteins
#for ( i in 1:650 ){
#  points(i,rawdata_sorted[i,"hsv_avgzscore"], pch=20, cex=1.5, col="olivedrab")
#}
#
#add known substrates
for ( i in 1:length(ICP0subs) ){
  currsub<-ICP0subs[i]
  points(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"hsv_avgzscore"], pch=20, cex=1.25, col="red")
  #text(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"hsv_avgzscore"], labels=rawdata_sorted[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
}
#
#add zinc finger proteins
#for ( i in 1:length(ZincFingerProts) ){
#  currsub<-ZincFingerProts[i]
#  points(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"hsv_avgzscore"], pch=20, cex=1.25, col="red")
#  text(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"hsv_avgzscore"], labels=rawdata_sorted[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
#}
#
axis(1,c(seq(minprotnum,maxprotnum,by=100)),col="black",cex.axis=1,lwd=3)
mtext("proteins (ordered by zscore)",1,line=2.5)
axis(2,c(seq(minzscore,maxzscore,by=0.5)),col="black",cex.axis=1,lwd=3)
mtext("mock zscore",2,line=2.5)
dev.off()
#
pdf(file.path(".", analysisdir, "dICP0_RawData_SnakePlot.pdf"), height=6, width=8)
rawdata_sorted<-rawdata[order(rawdata$dICP0_avgzscore, decreasing=TRUE), ]
nazscores <- apply(rawdata_sorted, 1, function(x) is.na(x[["dICP0_avgzscore"]]))
rawdata_sorted<-rawdata_sorted[ !nazscores, ]
rawdata_sorted$rownum<-1:nrow(rawdata_sorted)
minzscore<-floor(min(rawdata_sorted$dICP0_avgzscore, na.rm=TRUE))
maxzscore<-ceiling(max(rawdata_sorted$dICP0_avgzscore, na.rm=TRUE))
#minprotnum<-floor(min(rawdata_sorted$rownum, na.rm=TRUE))
minprotnum<-0
maxprotnum<-ceiling(max(rawdata_sorted$rownum, na.rm=TRUE))
par(bty='l',lwd=3)
plot(rawdata_sorted$dICP0_avgzscore,
     pch=20,
     cex=1.5,
     col="grey50",
     ylim=c(minzscore,maxzscore),
     xlab="", ylab="", xaxt='n', yaxt='n'
)
#color first 650 high abundance proteins
#for ( i in 1:650 ){
#  points(i,rawdata_sorted[i,"dICP0_avgzscore"], pch=20, cex=1.5, col="salmon")
#}
#
#add known substrates
for ( i in 1:length(ICP0subs) ){
  currsub<-ICP0subs[i]
  points(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"dICP0_avgzscore"], pch=20, cex=1.25, col="red")
  #text(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"dICP0_avgzscore"], labels=rawdata_sorted[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
}
#
#add zinc finger proteins
#for ( i in 1:length(ZincFingerProts) ){
#  currsub<-ZincFingerProts[i]
#  points(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"dICP0_avgzscore"], pch=20, cex=1.25, col="red")
#  text(rawdata_sorted[currsub,"rownum"],rawdata_sorted[currsub,"dICP0_avgzscore"], labels=rawdata_sorted[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
#}
#
axis(1,c(seq(minprotnum,maxprotnum,by=100)),col="black",cex.axis=1,lwd=3)
mtext("proteins (ordered by zscore)",1,line=2.5)
axis(2,c(seq(minzscore,maxzscore,by=0.5)),col="black",cex.axis=1,lwd=3)
mtext("mock zscore",2,line=2.5)
dev.off()
#
##########raw data snake plots


#####rawdata_TMP <- data.frame(rawdata)
#####rawdata_TMP[, "mock_avgzscore"][is.na(rawdata_TMP[, "mock_avgzscore"])] <- -4
#####rawdata_TMP[, "hsv_avgzscore"][is.na(rawdata_TMP[, "hsv_avgzscore"])] <- -4
#####rawdata_TMP[, "dICP0_avgzscore"][is.na(rawdata_TMP[, "dICP0_avgzscore"])] <- -4
#####
######raw: print all data
#####outfile=file.path(".", analysisdir, "HSVdICP0iPOND_RawDataAnalysis_AllDataWithImputedZscores_20200817.txt")
#####write.table(rawdata_TMP, file=outfile, sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
######


pdf(file.path(".", analysisdir, "HSVwtHSVdICP0_RawData_zscores.pdf"))
par(bty='l',lwd=3)
xmin <- floor(min(as.numeric(rawdata$hsv_modifiedzscore), na.rm=TRUE))
xmax <- ceiling(max(as.numeric(rawdata$hsv_modifiedzscore), na.rm=TRUE))
ymin <- floor(min(as.numeric(rawdata$dICP0_modifiedzscore), na.rm=TRUE))
ymax <- ceiling(max(as.numeric(rawdata$dICP0_modifiedzscore), na.rm=TRUE))
#xmin <- -2
#xmax <- 0
#ymin <- -4.5
#ymax <- -3.5
plot(rawdata$hsv_modifiedzscore, rawdata$dICP0_modifiedzscore,    
     pch=20,
     cex=1.2,
     col="grey70",
     xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     xlab="",ylab="",xaxt='n',yaxt='n'
)
#text(rawdata$hsv_modifiedzscore,rawdata$dICP0_modifiedzscore, labels=rawdata$geneID, cex=0.5, pos=4, col="red", font=2)
#
#add known substrates
for ( i in 1:length(ICP0subs) ){
  currsub<-ICP0subs[i]
  points(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], pch=20, cex=1.5, col="red")
  #text(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], labels=rawdata[currsub,"geneID"], cex=0.5, pos=2, col="red", font=2)
}
#
#add zinc finger proteins
for ( i in 1:length(ZincFingerProts) ){
  currsub<-ZincFingerProts[i]
  points(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], pch=20, cex=1.5, col="purple")
  #text(rawdata[currsub,"hsv_modifiedzscore"],rawdata[currsub,"dICP0_modifiedzscore"], labels=rawdata[currsub,"geneID"], cex=0.5, pos=1, col="purple", font=2)
}
#
axis(1,c(seq(xmin,xmax,by=1)),col="black",cex.axis=1,lwd=3)
mtext("HSV wt z-score",1,line=2.5)
axis(2,c(seq(ymin,ymax,by=0.5)),col="black",cex.axis=1,lwd=3)
mtext("HSV dICP0 z-score",2,line=2.5)
dev.off()


##########raw data volcano plots
#
##mock vs HSV WT volcano plot
pdf(file.path(".", analysisdir, "MockHSV_RawData_HSVWTdICP0enrichment_VolcanoPlot.pdf"))
minfc <- floor(min(as.numeric(rawdata$mock_hsv_avgfc), na.rm=TRUE))
maxfc <- ceiling(max(as.numeric(rawdata$mock_hsv_avgfc), na.rm=TRUE))
minpval <- floor(min(as.numeric(rawdata$mock_hsv_neglog10pval), na.rm=TRUE))
maxpval <- ceiling(max(as.numeric(rawdata$mock_hsv_neglog10pval), na.rm=TRUE))
#minfc <- -6
#maxfc <- 6
#minpval <- 0
#maxpval <- 3.5
par(bty='l',lwd=3)
plot(rawdata$mock_hsv_avgfc, rawdata$mock_hsv_neglog10pval,    
     pch=20,
     cex=1.2,
     #col="grey50",
     col=ifelse(rawdata$mock_hsv_neglog10pval<1.3, "grey75", "grey50"),
     xlim=c(minfc,maxfc), ylim=c(minpval,maxpval),
     xlab="",ylab="",xaxt='n',yaxt='n'
)
#
#add known substrates
#for ( i in 1:length(ICP0subs) ){
#  currsub<-ICP0subs[i]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.35, col="red")
#  text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
#}
#
#add proteins enriched on HSV WT compared to mock
#for ( i in 1:length(hsv_enriched_mock_prots_raw) ){
#  currsub<-hsv_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.2, col="purple")
#  text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="purple", font=2)
#}
#
#add proteins enriched on dICP0 compared to mock
#for ( i in 1:length(dICP0_enriched_mock_prots_raw) ){
#  currsub<-dICP0_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.2, col="olivedrab")
#  text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#}
#
#add proteins enriched on hsv compared to mock and dICP0 compared to mock
#for ( i in 1:length(hsv_dICP0_enriched_mock_prots_raw) ){
#  currsub<-hsv_dICP0_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.2, col="skyblue")
#  text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="skyblue", font=2)
#}
#
#add proteins with high abundance on HSV WT and dICP0
#for ( i in 1:length(rawdata_highabundance_HSVdCIP0intersect) ){
#  currsub<-rawdata_highabundance_HSVdCIP0intersect[i]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.2, col="skyblue")
#  #text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="skyblue", font=2)
#}
#
#add proteins with high abundance on HSV dICP0 but not WT
#for ( i in 1:length(rawdata_highabundance_ICP0sub_inhibit) ){
#  currsub<-rawdata_highabundance_ICP0sub_inhibit[i]
#  currfc<-rawdata[currsub,"mock_hsv_avgfc"]
#  currpval<-rawdata[currsub,"mock_hsv_neglog10pval"]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.35, col="red")
#  if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
#    text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
#  }
#}
#
#add proteins with high abundance on HSV WT but not dICP0
#for ( i in 1:length(rawdata_highabundance_ICP0sub_recruit) ){
#  currsub<-rawdata_highabundance_ICP0sub_recruit[i]
#  currfc<-rawdata[currsub,"mock_hsv_avgfc"]
#  currpval<-rawdata[currsub,"mock_hsv_neglog10pval"]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.35, col="olivedrab")
#  if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
#    text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#  }
#}
#
#add proteins enriched on HSV WT compared to dICP0
for ( i in 1:length(hsv_enriched_dICP0_prots_raw) ){
  currsub<-hsv_enriched_dICP0_prots_raw[i]
  currfc<-rawdata[currsub,"mock_hsv_avgfc"]
  currpval<-rawdata[currsub,"mock_hsv_neglog10pval"]
  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.35, col="gold")
  #if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
  #  text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="gold", font=2)
  #}
}
#
#add proteins enriched on HSV dICP0 compared to WT
for ( i in 1:length(dICP0_enriched_hsv_prots_raw) ){
  currsub<-dICP0_enriched_hsv_prots_raw[i]
  currfc<-rawdata[currsub,"mock_hsv_avgfc"]
  currpval<-rawdata[currsub,"mock_hsv_neglog10pval"]
  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.35, col="dodgerblue")
  #if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
  #  text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
  #}
}
#
#add clustered proteins
#for ( i in 1:length(PCAsubClusUnion) ){
#  currsub<-PCAsubClusUnion[i]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.2, col="dodgerblue")
#  #  text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
#}
#
#add clustered proteins
#for ( i in 1:length(PCAsubClusIntersect) ){
#  currsub<-PCAsubClusIntersect[i]
#  points(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.2, col="navy")
#  #text(rawdata[currsub,"mock_hsv_avgfc"],rawdata[currsub,"mock_hsv_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="navy", font=2)
#}
#add K11 subs
#for ( i in 1:length(K11subs) ){
#  currsub<-K11subs[i]
#  if ( currsub %in% rawdata$geneID  ){
#    points(rawdata[rawdata$geneID==currsub,"mock_hsv_avgfc"],rawdata[rawdata$geneID==currsub,"mock_hsv_neglog10pval"], pch=20, cex=1.25, col="olivedrab")
#    text(rawdata[rawdata$geneID==currsub,"mock_hsv_avgfc"],rawdata[rawdata$geneID==currsub,"mock_hsv_neglog10pval"], labels=rawdata[rawdata$geneID==currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#  }
#}
#
axis(1,c(seq(minfc,maxfc,by=1)),col="black",cex.axis=1,lwd=3)
mtext("HSV-mock log2 fold change",1,line=2.5)
axis(2,c(seq(minpval,maxpval,by=0.5)),col="black",cex.axis=1,lwd=3)
mtext("HSV-mock -log10 pval",2,line=2.5)
dev.off()
#
##mock vs dICP0 volcano plot
pdf(file.path(".", analysisdir, "MockdICP0_RawData_HSVWTdICP0enrichment_VolcanoPlot.pdf"))
minfc <- floor(min(as.numeric(rawdata$mock_dICP0_avgfc), na.rm=TRUE))
maxfc <- ceiling(max(as.numeric(rawdata$mock_dICP0_avgfc), na.rm=TRUE))
minpval <- floor(min(as.numeric(rawdata$mock_dICP0_neglog10pval), na.rm=TRUE))
maxpval <- ceiling(max(as.numeric(rawdata$mock_dICP0_neglog10pval), na.rm=TRUE))
#minfc <- -5
#maxfc <- 5
#minpval <- 0
#maxpval <- 3
par(bty='l',lwd=3)
plot(rawdata$mock_dICP0_avgfc, rawdata$mock_dICP0_neglog10pval,
     pch=20,
     cex=1.2,
     #col="grey50",
     col=ifelse(rawdata$mock_dICP0_neglog10pval<1.3, "grey75", "grey50"),
     xlim=c(minfc,maxfc), ylim=c(minpval,maxpval),
     xlab="",ylab="",xaxt='n',yaxt='n'
)
#
#add known substrates
#for ( i in 1:length(ICP0subs) ){
#  currsum<-ICP0subs[i]
#  points(rawdata[currsum,"mock_dICP0_avgfc"],rawdata[currsum,"mock_dICP0_neglog10pval"], pch=20, cex=1.35, col="red")
#  text(rawdata[currsum,"mock_dICP0_avgfc"],rawdata[currsum,"mock_dICP0_neglog10pval"], labels=rawdata[currsum,"geneID"], cex=0.75, pos=2, col="red", font=2)
#}
#
#add proteins enriched on HSV WT compared to mock
#for ( i in 1:length(hsv_enriched_mock_prots_raw) ){
#  currsub<-hsv_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.2, col="purple")
#  text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="purple", font=2)
#}
#
#add proteins enriched on dICP0 compared to mock
#for ( i in 1:length(dICP0_enriched_mock_prots_raw) ){
#  currsub<-dICP0_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.2, col="olivedrab")
#  text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#}
#
#add proteins enriched on hsv compared to mock and dICP0 compared to mock
#for ( i in 1:length(hsv_dICP0_enriched_mock_prots_raw) ){
#  currsub<-hsv_dICP0_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.2, col="skyblue")
#  text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="skyblue", font=2)
#}
#
#add proteins with high abundance on HSV WT and dICP0
#for ( i in 1:length(rawdata_highabundance_HSVdCIP0intersect) ){
#  currsub<-rawdata_highabundance_HSVdCIP0intersect[i]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.2, col="skyblue")
#  #text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="skyblue", font=2)
#}
#
#add proteins with high abundance on HSV dICP0 but not WT
#for ( i in 1:length(rawdata_highabundance_ICP0sub_inhibit) ){
#  currsub<-rawdata_highabundance_ICP0sub_inhibit[i]
#  currfc<-rawdata[currsub,"mock_dICP0_avgfc"]
#  currpval<-rawdata[currsub,"mock_dICP0_neglog10pval"]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.35, col="red")
#  if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
#    text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
#  }
#}
#
#add proteins with high abundance on HSV WT but not dICP0
#for ( i in 1:length(rawdata_highabundance_ICP0sub_recruit) ){
#  currsub<-rawdata_highabundance_ICP0sub_recruit[i]
#  currfc<-rawdata[currsub,"mock_dICP0_avgfc"]
#  currpval<-rawdata[currsub,"mock_dICP0_neglog10pval"]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.35, col="olivedrab")
#  if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
#    text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#  }
#}
#
#add proteins enriched on HSV WT compared to dICP0
for ( i in 1:length(hsv_enriched_dICP0_prots_raw) ){
  currsub<-hsv_enriched_dICP0_prots_raw[i]
  currfc<-rawdata[currsub,"mock_dICP0_avgfc"]
  currpval<-rawdata[currsub,"mock_dICP0_neglog10pval"]
  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.35, col="gold")
  #if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
  #  text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="gold", font=2)
  #}
}
#
#add proteins enriched on HSV dICP0 compared to WT
for ( i in 1:length(dICP0_enriched_hsv_prots_raw) ){
  currsub<-dICP0_enriched_hsv_prots_raw[i]
  currfc<-rawdata[currsub,"mock_dICP0_avgfc"]
  currpval<-rawdata[currsub,"mock_dICP0_neglog10pval"]
  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.35, col="dodgerblue")
  #if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
  #  text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
  #}
}
#
#
#add clustered proteins
#for ( i in 1:length(PCAsubClusUnion) ){
#  currsub<-PCAsubClusUnion[i]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.25, col="dodgerblue")
#  #  text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
#}
#
#add clustered proteins
#for ( i in 1:length(PCAsubClusIntersect) ){
#  currsub<-PCAsubClusIntersect[i]
#  points(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.25, col="navy")
#  text(rawdata[currsub,"mock_dICP0_avgfc"],rawdata[currsub,"mock_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="navy", font=2)
#}
#add K11 subs
#for ( i in 1:length(K11subs) ){
#  currsub<-K11subs[i]
#  if ( currsub %in% rawdata$geneID  ){
#    points(rawdata[rawdata$geneID==currsub,"mock_dICP0_avgfc"],rawdata[rawdata$geneID==currsub,"mock_dICP0_neglog10pval"], pch=20, cex=1.25, col="olivedrab")
#    text(rawdata[rawdata$geneID==currsub,"mock_dICP0_avgfc"],rawdata[rawdata$geneID==currsub,"mock_dICP0_neglog10pval"], labels=rawdata[rawdata$geneID==currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#  }
#}
#
axis(1,c(seq(minfc,maxfc,by=1)),col="black",cex.axis=1,lwd=3)
mtext("dICP0-mock log2 fold change",1,line=2.5)
axis(2,c(seq(minpval,maxpval,by=0.5)),col="black",cex.axis=1,lwd=3)
mtext("dICP0-mock -log10 pval",2,line=2.5)
dev.off()
#
##HSV WT vs dICP0 volcano plot
pdf(file.path(".", analysisdir, "HSVdICP0_RawData_HSVWTdICP0enrichment_nolabel_VolcanoPlot.pdf"))
minfc <- floor(min(as.numeric(rawdata$hsv_dICP0_avgfc), na.rm=TRUE))
maxfc <- ceiling(max(as.numeric(rawdata$hsv_dICP0_avgfc), na.rm=TRUE))
minpval <- floor(min(as.numeric(rawdata$hsv_dICP0_neglog10pval), na.rm=TRUE))
maxpval <- ceiling(max(as.numeric(rawdata$hsv_dICP0_neglog10pval), na.rm=TRUE))
#minfc <- -6
#maxfc <- 6
#minpval <- 0
#maxpval <- 3
par(bty='l',lwd=3)
plot(rawdata$hsv_dICP0_avgfc, rawdata$hsv_dICP0_neglog10pval,
     pch=20,
     cex=1.2,
     #col="grey50",
     col=ifelse(rawdata$hsv_dICP0_neglog10pval<1.3, "grey75", "grey50"),
     xlim=c(minfc,maxfc), ylim=c(minpval,maxpval),
     xlab="",ylab="",xaxt='n',yaxt='n'
)
#
#add known substrates
#for ( i in 1:length(ICP0subs) ){
#  currsum<-ICP0subs[i]
#  points(rawdata[currsum,"hsv_dICP0_avgfc"],rawdata[currsum,"hsv_dICP0_neglog10pval"], pch=20, cex=1.35, col="red")
#  text(rawdata[currsum,"hsv_dICP0_avgfc"],rawdata[currsum,"hsv_dICP0_neglog10pval"], labels=rawdata[currsum,"geneID"], cex=0.75, pos=2, col="red", font=2)
#}
#
#add proteins enriched on HSV WT compared to mock
#for ( i in 1:length(hsv_enriched_mock_prots_raw) ){
#  currsub<-hsv_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.2, col="purple")
#  text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="purple", font=2)
#}
#
#add proteins enriched on dICP0 compared to mock
#for ( i in 1:length(dICP0_enriched_mock_prots_raw) ){
#  currsub<-dICP0_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.2, col="olivedrab")
#  text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#}
#
#add proteins enriched on hsv compared to mock and dICP0 compared to mock
#for ( i in 1:length(hsv_dICP0_enriched_mock_prots_raw) ){
#  currsub<-hsv_dICP0_enriched_mock_prots_raw[i]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.2, col="skyblue")
#  text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="skyblue", font=2)
#}
#
#add proteins with high abundance on HSV WT and dICP0
#for ( i in 1:length(rawdata_highabundance_HSVdCIP0intersect) ){
#  currsub<-rawdata_highabundance_HSVdCIP0intersect[i]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.2, col="skyblue")
#  #text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="skyblue", font=2)
#}
#
#add proteins with high abundance on HSV dICP0 but not WT
#for ( i in 1:length(rawdata_highabundance_ICP0sub_inhibit) ){
#  currsub<-rawdata_highabundance_ICP0sub_inhibit[i]
#  currfc<-rawdata[currsub,"hsv_dICP0_avgfc"]
#  currpval<-rawdata[currsub,"hsv_dICP0_neglog10pval"]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.35, col="red")
#  if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
#    text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="red", font=2)
#  }
#}
#
#add proteins with high abundance on HSV WT but not dICP0
#for ( i in 1:length(rawdata_highabundance_ICP0sub_recruit) ){
#  currsub<-rawdata_highabundance_ICP0sub_recruit[i]
#  currfc<-rawdata[currsub,"hsv_dICP0_avgfc"]
#  currpval<-rawdata[currsub,"hsv_dICP0_neglog10pval"]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.35, col="olivedrab")
#  if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
#    text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#  }
#}
#
#add proteins enriched on HSV WT compared to dICP0
for ( i in 1:length(hsv_enriched_dICP0_prots_raw) ){
  currsub<-hsv_enriched_dICP0_prots_raw[i]
  currfc<-rawdata[currsub,"hsv_dICP0_avgfc"]
  currpval<-rawdata[currsub,"hsv_dICP0_neglog10pval"]
  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.35, col="gold")
  #text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="gold", font=2)
  #if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
  #  text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="gold", font=2)
  #}
}
#
#add proteins enriched on HSV dICP0 compared to WT
for ( i in 1:length(dICP0_enriched_hsv_prots_raw) ){
  currsub<-dICP0_enriched_hsv_prots_raw[i]
  currfc<-rawdata[currsub,"hsv_dICP0_avgfc"]
  currpval<-rawdata[currsub,"hsv_dICP0_neglog10pval"]
  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.35, col="dodgerblue")
  #text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
  #if ( (!is.na(currfc) & (currfc > 0.5 | currfc < -0.5)) & (!is.na(currpval) & (currpval > 1)) ){
  #  text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
  #}
}
#
#add clustered proteins
#for ( i in 1:length(PCAsubClusUnion) ){
#  currsub<-PCAsubClusUnion[i]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.25, col="dodgerblue")
#  #  text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="dodgerblue", font=2)
#}
#
#add clustered proteins
#for ( i in 1:length(PCAsubClusIntersect) ){
#  currsub<-PCAsubClusIntersect[i]
#  points(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.25, col="navy")
#  text(rawdata[currsub,"hsv_dICP0_avgfc"],rawdata[currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[currsub,"geneID"], cex=0.75, pos=2, col="navy", font=2)
#}
#add K11 subs
#for ( i in 1:length(K11subs) ){
#  currsub<-K11subs[i]
#  if ( currsub %in% rawdata$geneID  ){
#    points(rawdata[rawdata$geneID==currsub,"hsv_dICP0_avgfc"],rawdata[rawdata$geneID==currsub,"hsv_dICP0_neglog10pval"], pch=20, cex=1.25, col="olivedrab")
#    text(rawdata[rawdata$geneID==currsub,"hsv_dICP0_avgfc"],rawdata[rawdata$geneID==currsub,"hsv_dICP0_neglog10pval"], labels=rawdata[rawdata$geneID==currsub,"geneID"], cex=0.75, pos=2, col="olivedrab", font=2)
#  }
#}
#
axis(1,c(seq(minfc,maxfc,by=1)),col="black",cex.axis=1,lwd=3)
mtext("dICP0-HSV log2 fold change",1,line=2.5)
axis(2,c(seq(minpval,maxpval,by=0.5)),col="black",cex.axis=1,lwd=3)
mtext("dICP0-HSV -log10 pval",2,line=2.5)
dev.off()
#
##########end raw data volcano plots


##########heatmap of mock, HSV, and dICP0 raw z-scores for high abundance proteins on HSV WT or dICP0
#
MORC3network<-c("Q14149", "P10155", "Q14839", "Q71DI3", "P23246", "Q12888", "Q9UBE0", "Q52LJ0", "Q01780", "Q8N1G2", "Q9UPN9", "Q9Y3I0", "P52272", "O94822")
NACC1network<-c("O00505", "O75182", "P61513", "Q96RE7", "Q96ST3", "Q8IU81", "O14497", "Q13618", "P04632", "Q99459", "O00629", "Q13148", "P35268", "P15880")
ZBTB20network<-c("P12109", "Q9HC78", "O43670")
ZEB1network<-c("P09038", "P51532", "P19022", "Q9HAV4", "P08670", "P28482", "P08134", "P09651", "P26358", "Q13315", "P40763", "P37275", "Q13363", "Q15436", "P02751", "Q13547", "P08758", "P56545", "P00352", "Q93009", "P15407", "Q12888", "P63165", "O60341", "Q07157", "P26447", "P04406", "P35222", "P16070")
ZNF644network<-c("Q14966", "Q9Y232", "O15460", "Q9H582", "Q96KQ7", "Q8IVL6", "Q9Y376", "O95218", "O95785")
#rawdata_HSV_high <- rawdata[row.names(rawdata)%in%rawdata_highabundance_ICP0sub_recruit,]
rawdata_HSV_high <- rawdata[row.names(rawdata)%in%ZEB1network,]
rawdata_HSV_high <- rawdata_HSV_high[order(rawdata_HSV_high$dICP0_avgzscore, decreasing=TRUE),]
rawdata_dICP0_high <- rawdata[row.names(rawdata)%in%rawdata_highabundance_ICP0sub_inhibit,]
rawdata_dICP0_high <- rawdata_dICP0_high[order(rawdata_dICP0_high$dICP0_avgzscore, decreasing=TRUE),]
#rawdata_enriched <- Reduce(rbind, list(rawdata_HSV_high, rawdata_dICP0_high))
#rawdata_enriched <- rawdata_dICP0_high
rawdata_enriched <- rawdata_HSV_high
dim(rawdata_HSV_high)
dim(rawdata_dICP0_high)
dim(rawdata_enriched)
#
prots<-row.names(rawdata_enriched)
genes<-rawdata_enriched$geneID
print(genes)
print(prots)
#
#rawdata_enriched <- rawdata_enriched[order(rawdata_enriched$hsv_avgzscore, decreasing=TRUE),]
#
rawdata_enriched_zscores <- rawdata_enriched[,c("hsv_avgzscore", "dICP0_avgzscore")]
#
rawdata_enriched_zscores_mat <- as.matrix(rawdata_enriched_zscores)
rawdata_enriched_zscores_mat <- apply(rawdata_enriched_zscores_mat,2,as.numeric)
rawdata_enriched_zscores_mat<- round(rawdata_enriched_zscores_mat, 4)
row.names(rawdata_enriched_zscores_mat) <- rawdata_enriched$geneID
head(rawdata_enriched_zscores_mat)
#
##heatmap of avg intensities
pdf(file.path(".", analysisdir, ("MockHSVdICP0_avgZscore_ZEB1interactome_heatmap.pdf")))
col_palette <- colorRampPalette(c("yellow", "orange", "red"))(n = 98)
#
#hmcol<-brewer.pal(9,"YlGnBu")
hmcol<-brewer.pal(9,"YlOrRd")
#hmcol<-brewer.pal(9, "RdYlBu")
#hmcol<-brewer.pal(9, "RdBu")
#hmcol<-brewer.pal(9, "Reds")
heatmap.2(rawdata_enriched_zscores_mat,
          notecex=1.5,
          notecol="grey80",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #dendrogram="row",    # only draw a row dendrogram
          #Rowv=TRUE,
          dendrogram="none",
          Rowv="NA",
          Colv="NA",            # turn off column clustering
          #col=bluered(100),     #col_palette or col=bluered(100),
          #col=col_palette,
          #breaks=col_breaks,
          col=hmcol,
          na.color = "grey",
          scale="none",
          #labRow=CullinData$GeneId.x,
          labCol=c("hsv", "dICP0"),
          margins=c(6,20),
          cexRow=0.5,
          cexCol=1,
          vline=NULL,
          hline=NULL,
          key=TRUE,
          keysize=1.0,
          sepwidth=c(0.0, 0.01),  # width of the borders
          sepcolor='white',
          colsep=seq(0,ncol(rawdata_enriched_zscores_mat),by=1),
          #rowsep=seq(0,nrow(rawdata_enriched_zscores_mat),by=1)
)
dev.off()
#
##########end heatmap of mock, HSV, and dICP0 raw z-scores for high abundance proteins on HSV WT or dICP0





########## this isn't used ## make heatmap of mock, HSV, and dICP0 z-scores for high abundance proteins on HSV WT or dICP0
rawdata_highabundance_prot<-union(rawdata_HSV_highabundance_prot, rawdata_dICP0_highabundance_prot)
#rawdata_HSVenriched <- rawdata[row.names(rawdata)%in%rawdata_HSV_highabundance_prot,]
#rawdata_HSVenriched <- rawdata_HSVenriched[order(rawdata_HSVenriched$hsv_avgzscore, decreasing=TRUE),]
#dim(rawdata_HSVenriched)
#rawdata_dICP0enriched <- rawdata[row.names(rawdata)%in%rawdata_dICP0_highabundance_prot,]
#rawdata_dICP0enriched <- rawdata_dICP0enriched[order(rawdata_dICP0enriched$dICP0_avgzscore, decreasing=TRUE),]
#dim(rawdata_dICP0enriched)
#rawdata_enriched <- Reduce(rbind, list(rawdata_HSVenriched, rawdata_dICP0enriched))
rawdata_enriched <- rawdata[row.names(rawdata)%in%rawdata_highabundance_prot,]
dim(rawdata_enriched)
#
prots<-row.names(rawdata_enriched)
genes<-rawdata_enriched$geneID
#
#rawdata_enriched <- rawdata_enriched[order(rawdata_enriched$hsv_avgzscore, decreasing=TRUE),]
#
narows <- apply(rawdata_enriched, 1, function(x) ifelse((is.na(x["hsv_avgzscore"]) | is.na(x["dICP0_avgzscore"])),TRUE,FALSE))
rawdata_enriched <- rawdata_enriched[ !narows, ]
rawdata_enriched_zscores <- rawdata_enriched[,c("hsv_avgzscore", "dICP0_avgzscore")]
dim(rawdata_enriched_zscores)
#narows <- apply(rawdata_enriched_zscores, 1, function(x) ifelse((is.na(x["mock_avgzscore"]) | is.na(x["hsv_avgzscore"]) | is.na(x["dICP0_avgzscore"])),TRUE,FALSE))
#rawdata_enriched_zscores <- rawdata_enriched_zscores[ !narows, ]
dim(rawdata_enriched_zscores)
#
rawdata_enriched_zscores_mat <- as.matrix(rawdata_enriched_zscores)
head(rawdata_enriched_zscores_mat)
rawdata_enriched_zscores_mat <- apply(rawdata_enriched_zscores_mat,2,as.numeric)
rawdata_enriched_zscores_mat<- round(rawdata_enriched_zscores_mat, 2)
row.names(rawdata_enriched_zscores_mat) <- rawdata_enriched$geneID
head(rawdata_enriched_zscores_mat)
#
##heatmap of avg intensities
pdf(file.path(".", analysisdir, ("MockHSVdICP0_HSV-dICP0_highabundanceprots_heatmap.pdf")))
col_palette <- colorRampPalette(c("yellow", "green", "blue"))(n = 98)
#
hmcol<-brewer.pal(9,"YlGnBu")
#hmcol<-brewer.pal(9,"YlOrRd")
#hmcol<-brewer.pal(9, "RdYlBu")
#hmcol<-brewer.pal(9, "RdBu")
#hmcol<-brewer.pal(9, "Reds")
heatmap.2(rawdata_enriched_zscores_mat,
          notecex=1.5,
          notecol="grey80",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="row",    # only draw a row dendrogram
          Rowv=TRUE,
          #dendrogram="none",
          #Rowv="NA",
          Colv="NA",            # turn off column clustering
          #col=bluered(100),     #col_palette or col=bluered(100),
          #col=col_palette,
          #breaks=col_breaks,
          col=hmcol,
          na.color = "grey",
          scale="none",
          #labRow=CullinData$GeneId.x,
          labCol=c("mock", "hsv", "dICP0"),
          margins=c(6,20),
          cexRow=0.2,
          cexCol=1,
          vline=NULL,
          hline=NULL,
          key=TRUE,
          keysize=1.0,
          sepwidth=c(0.0, 0.01),  # width of the borders
          sepcolor='white',
          colsep=seq(0,ncol(rawdata_enriched_zscores_mat),by=1)
          #rowsep=seq(0,nrow(rawdata_enriched_zscores_mat),by=1)
)
dev.off()
#
########## this isn't used ## make heatmap of mock, HSV, and dICP0 z-scores for high abundance proteins on HSV WT or dICP0












#make heatmap of mock, HSV, and dICP0 z-scores for proteins enriched on HSV WT compared to mock
rawdata_HSVdICP0_enriched <- rawdata[row.names(rawdata)%in%hsv_dICP0_enriched_mock_prots_raw,]
rawdata_HSVdICP0_enriched <- rawdata_HSVdICP0_enriched[order(rawdata_HSVdICP0_enriched$hsv_avgzscore, decreasing=TRUE),]
rawdata_HSVenriched <- rawdata[row.names(rawdata)%in%hsv_enriched_dICP0_notenriched_mock_prots_raw,]
rawdata_HSVenriched <- rawdata_HSVenriched[order(rawdata_HSVenriched$hsv_avgzscore, decreasing=TRUE),]
rawdata_dICP0enriched <- rawdata[row.names(rawdata)%in%dICP0_enriched_hsv_notenriched_mock_prots_raw,]
rawdata_dICP0enriched <- rawdata_dICP0enriched[order(rawdata_dICP0enriched$dICP0_avgzscore, decreasing=TRUE),]
rawdata_enriched <- Reduce(rbind, list(rawdata_HSVdICP0_enriched, rawdata_HSVenriched, rawdata_dICP0enriched))
dim(rawdata_HSVdICP0_enriched)
dim(rawdata_HSVenriched)
dim(rawdata_dICP0enriched)
dim(rawdata_enriched)
#
prots<-row.names(rawdata_enriched)
genes<-rawdata_enriched$geneID
print(genes)
print(prots)
#
#rawdata_enriched <- rawdata_enriched[order(rawdata_enriched$hsv_avgzscore, decreasing=TRUE),]
#
rawdata_HSVenriched_zscores <- rawdata_enriched[,c("mock_avgzscore", "hsv_avgzscore", "dICP0_avgzscore")]
#
rawdata_enriched_zscores_mat <- as.matrix(rawdata_HSVenriched_zscores)
rawdata_enriched_zscores_mat <- apply(rawdata_enriched_zscores_mat,2,as.numeric)
rawdata_enriched_zscores_mat<- round(rawdata_enriched_zscores_mat, 2)
row.names(rawdata_enriched_zscores_mat) <- rawdata_enriched$geneID
head(rawdata_enriched_zscores_mat)
#
##heatmap of avg intensities
pdf(file.path(".", analysisdir, ("MockHSVdICP0_avgZscore_HSVenrichedMock-dICP0enrichedMock_heatmap.pdf")))
col_palette <- colorRampPalette(c("yellow", "green", "blue"))(n = 98)
#
hmcol<-brewer.pal(9,"YlGnBu")
#hmcol<-brewer.pal(9,"YlOrRd")
#hmcol<-brewer.pal(9, "RdYlBu")
#hmcol<-brewer.pal(9, "RdBu")
#hmcol<-brewer.pal(9, "Reds")
heatmap.2(rawdata_enriched_zscores_mat,
          notecex=1.5,
          notecol="grey80",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #dendrogram="row",    # only draw a row dendrogram
          #Rowv=TRUE,
          dendrogram="none",
          Rowv="NA",
          Colv="NA",            # turn off column clustering
          #col=bluered(100),     #col_palette or col=bluered(100),
          #col=col_palette,
          #breaks=col_breaks,
          col=hmcol,
          na.color = "grey",
          scale="none",
          #labRow=CullinData$GeneId.x,
          labCol=c("mock", "hsv", "dICP0"),
          margins=c(6,20),
          cexRow=0.6,
          cexCol=1,
          vline=NULL,
          hline=NULL,
          key=TRUE,
          keysize=1.0,
          sepwidth=c(0.0, 0.01),  # width of the borders
          sepcolor='white',
          colsep=seq(0,ncol(rawdata_enriched_zscores_mat),by=1),
          #colsep=0:ncol(rawdata_enriched_zscores_mat),
          rowsep=seq(0,nrow(rawdata_enriched_zscores_mat),by=1)
)
dev.off()






##volcano plots
#
plotvolcanoplot(imputeddata)
#
##end volcano plots






#make heatmap of mock, HSV, and dICP0 z-scores
dim(imputeddata)
allgenes<-vector()
allgenes<-imputeddata$geneID
allgenes<-as.character(allgenes)
write(allgenes, "alliPONDgenes.txt", ncolumns=1)
#imputeddata_sigenriched <- subset(imputeddata, (imputeddata$mock_hsv_TtestPval<0.05 | imputeddata$mock_dICP0_TtestPval<0.05 | imputeddata$hsv_dICP0_TtestPval<0.05))
#imputeddata_sigenriched <- subset(imputeddata, (imputeddata$hsv_dICP0_avgfc>0 & imputeddata$hsv_dICP0_TtestPval<0.05))
imputeddata_sigenriched <- subset(imputeddata, imputeddata$hsv_dICP0_avgfc>1)
head(imputeddata_sigenriched)
dim(imputeddata_sigenriched)
#
prots<-row.names(imputeddata_sigenriched)
genes<-imputeddata_sigenriched$geneID
print(genes)
print(prots)
#
##########imputeddata_sigenriched <- imputeddata_sigenriched[order(imputeddata_sigenriched$dICP0_avgzscore, decreasing=TRUE),]
#
imputeddata_zscores <- imputeddata_sigenriched[,c("mock_avgzscore", "hsv_avgzscore", "dICP0_avgzscore")]
#
imputeddata_zscores_mat <- as.matrix(imputeddata_zscores)
imputeddata_zscores_mat <- apply(imputeddata_zscores_mat,2,as.numeric)
imputeddata_zscores_mat<- round(imputeddata_zscores_mat, 2)
row.names(imputeddata_zscores_mat) <- imputeddata_sigenriched$geneID
head(imputeddata_zscores_mat)
#
##heatmap of avg intensities
print("heatmap")
pdf(file.path(".", analysisdir, ("MockHSVdICP0_avgZscore_HSVdICP0fcGT1_heatmap.pdf")))
col_palette <- colorRampPalette(c("yellow", "green", "blue"))(n = 98)
#
hmcol<-brewer.pal(9,"YlGnBu")
#hmcol<-brewer.pal(9,"YlOrRd")
#hmcol<-brewer.pal(9, "RdYlBu")
#hmcol<-brewer.pal(9, "RdBu")
#hmcol<-brewer.pal(9, "Reds")
heatmap.2(imputeddata_zscores_mat,
          #cellnote=KeGGdataRepIntensitiesMat,
          notecex=1.5,
          notecol="grey80",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="row",    # only draw a row dendrogram
          Rowv=TRUE,
          #dendrogram="none",
          #Rowv="NA",
          Colv="NA",            # turn off column clustering
          #col=bluered(100),     #col_palette or col=bluered(100),
          #col=col_palette,
          #breaks=col_breaks,
          col=hmcol,
          na.color = "grey",
          scale="none",
          #labRow=CullinData$GeneId.x,
          labCol=c("mock", "hsv", "dICP0"),
          margins=c(6,20),
          cexRow=0.2,
          cexCol=1,
          vline=NULL,
          hline=NULL,
          key=TRUE,
          keysize=1.0,
          sepwidth=c(0.0, 0.01),  # width of the borders
          sepcolor='white',
          colsep=seq(0,ncol(imputeddata_zscores_mat),by=1)
          #colsep=0:ncol(KeGGdataRepIntensitiesMat),
          #rowsep=0:nrow(imputeddata_zscores_mat)
)
dev.off()



stop("VOL")

#snake plot
snakeplot(imputeddata)



### redo these subroutines
##fold change vs abundance z-score
#plotfcabundance(repfc,imputeddata_zscores,expsigtest)
#

##individual replicate fold change distributions
#plotfoldchanges(mock_hsv_1_foldchanges, mock_hsv_2_foldchanges, mock_hsv_3_foldchanges, "MockHSV_FoldChangeDistribution.pdf")
#plotfoldchanges(mock_mut_1_foldchanges, mock_mut_2_foldchanges, mock_mut_3_foldchanges, "MockMut_FoldChangeDistribution.pdf")
#plotfoldchanges(hsv_mut_1_foldchanges, hsv_mut_2_foldchanges, hsv_mut_3_foldchanges, "HSVMut_FoldChangeDistribution.pdf")
#





stop("NN")



##PCA cluster analysis
#union
ALLclus<-vector()
ALLclus<-scan("IntersectClus.txt", what='complex')
print(ALLclus)

ICP0subs<-vector()
ICP0subs<-c("P78527", "P29590", "Q16666", "P46100")

###

ClustZscores<-cbind(imputeddata_zscores[ALLclus,"mock_avgzscore"], imputeddata_zscores[ALLclus,"hsv_avgzscore"], imputeddata_zscores[ALLclus,"dICP0_avgzscore"])
ClustZscores<-as.data.frame(ClustZscores, row.names=ALLclus)
colnames(ClustZscores)<-c("Mock.Zscore", "HSVwt.Zscore", "HSVdICP0.Zscore")
ClustZscores$ICP0sub<-


ClustZscoresLong<-melt(as.matrix(ClustZscores))
print(ClustZscores)
print(ClustZscoresLong)


linecol<-c("red","orange","yellow","green","blue","purple","grey","black","olivedrab","salmon","goldenrod","tomato")
ZscorePlot <- ggplot(ClustZscoresLong, aes(X2, value, group=factor(X1))) + geom_line(size = 1, colour="grey")

#ZscorePlot + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size=2))
ZscorePlot

stop("stopped")

if ( FALSE )
{
ClusData<-data.frame()
ALLclusdata<-imputeddata_zscores[row.names(imputeddata_zscores) %in% ALLclus, ]
ALLclusdata_MockMean<-mean(ALLclusdata$mock_avgzscore, na.rm=TRUE)
ALLclusdata_MockStdev<-sd(ALLclusdata$mock_avgzscore, na.rm=TRUE)
ALLclusdata_WtMean<-mean(ALLclusdata$hsv_avgzscore, na.rm=TRUE)
ALLclusdata_WtStdev<-sd(ALLclusdata$hsv_avgzscore, na.rm=TRUE)
ALLclusdata_MutMean<-mean(ALLclusdata$dICP0_avgzscore, na.rm=TRUE)
ALLclusdata_MutStdev<-sd(ALLclusdata$dICP0_avgzscore, na.rm=TRUE)
ClusData["ALLclus","MockMeanZscore"]<-ALLclusdata_MockMean
ClusData["ALLclus","MockStdevZscore"]<-ALLclusdata_MockStdev
ClusData["ALLclus","WtMeanZscore"]<-ALLclusdata_WtMean
ClusData["ALLclus","WtStdevZscore"]<-ALLclusdata_WtStdev
ClusData["ALLclus","MutMeanZscore"]<-ALLclusdata_MutMean
ClusData["ALLclus","MutStdevZscore"]<-ALLclusdata_MutStdev
ALLclusFCdata<-repfc[row.names(repfc) %in% ALLclus, ]
ALLclusFCdata_MockWtFcMean<-mean(ALLclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
ALLclusFCdata_MockWtFcStdev<-sd(ALLclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
ALLclusFCdata_MockMutFcMean<-mean(ALLclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
ALLclusFCdata_MockMutFcStdev<-sd(ALLclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
ALLclusFCdata_WtMutFcMean<-mean(ALLclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ALLclusFCdata_WtMutFcStdev<-sd(ALLclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ClusData["ALLclus","MockWtMeanFC"]<-ALLclusFCdata_MockWtFcMean
ClusData["ALLclus","MockWtStdevFC"]<-ALLclusFCdata_MockWtFcStdev
ClusData["ALLclus","MockMutMeanFC"]<-ALLclusFCdata_MockMutFcMean
ClusData["ALLclus","MockMutStdevFC"]<-ALLclusFCdata_MockMutFcStdev
ClusData["ALLclus","WtMutMeanFC"]<-ALLclusFCdata_WtMutFcMean
ClusData["ALLclus","WtMutStdevFC"]<-ALLclusFCdata_WtMutFcStdev
#USP7
USP7clus<-vector()
USP7clus<-scan("USP7clus.txt", what='complex')
USP7clusdata<-imputeddata_zscores[row.names(imputeddata_zscores) %in% USP7clus, ]
USP7clusdata_MockMean<-mean(USP7clusdata$mock_avgzscore, na.rm=TRUE)
USP7clusdata_MockStdev<-sd(USP7clusdata$mock_avgzscore, na.rm=TRUE)
USP7clusdata_WtMean<-mean(USP7clusdata$hsv_avgzscore, na.rm=TRUE)
USP7clusdata_WtStdev<-sd(USP7clusdata$hsv_avgzscore, na.rm=TRUE)
USP7clusdata_MutMean<-mean(USP7clusdata$dICP0_avgzscore, na.rm=TRUE)
USP7clusdata_MutStdev<-sd(USP7clusdata$dICP0_avgzscore, na.rm=TRUE)
ClusData["USP7clus","MockMeanZscore"]<-USP7clusdata_MockMean
ClusData["USP7clus","MockStdevZscore"]<-USP7clusdata_MockStdev
ClusData["USP7clus","WtMeanZscore"]<-USP7clusdata_WtMean
ClusData["USP7clus","WtStdevZscore"]<-USP7clusdata_WtStdev
ClusData["USP7clus","MutMeanZscore"]<-USP7clusdata_MutMean
ClusData["USP7clus","MutStdevZscore"]<-USP7clusdata_MutStdev
USP7clusFCdata<-repfc[row.names(repfc) %in% USP7clus, ]
USP7clusFCdata_MockWtFcMean<-mean(USP7clusFCdata$mock_hsv_avgfc, na.rm=TRUE)
USP7clusFCdata_MockWtFcStdev<-sd(USP7clusFCdata$mock_hsv_avgfc, na.rm=TRUE)
USP7clusFCdata_MockMutFcMean<-mean(USP7clusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
USP7clusFCdata_MockMutFcStdev<-sd(USP7clusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
USP7clusFCdata_WtMutFcMean<-mean(USP7clusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
USP7clusFCdata_WtMutFcStdev<-sd(USP7clusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ClusData["USP7clus","MockWtMeanFC"]<-USP7clusFCdata_MockWtFcMean
ClusData["USP7clus","MockWtStdevFC"]<-USP7clusFCdata_MockWtFcStdev
ClusData["USP7clus","MockMutMeanFC"]<-USP7clusFCdata_MockMutFcMean
ClusData["USP7clus","MockMutStdevFC"]<-USP7clusFCdata_MockMutFcStdev
ClusData["USP7clus","WtMutMeanFC"]<-USP7clusFCdata_WtMutFcMean
ClusData["USP7clus","WtMutStdevFC"]<-USP7clusFCdata_WtMutFcStdev
#PML
PMLclus<-vector()
PMLclus<-scan("PMLclus.txt", what='complex')
PMLclusdata<-imputeddata_zscores[row.names(imputeddata_zscores) %in% PMLclus, ]
PMLclusdata_MockMean<-mean(PMLclusdata$mock_avgzscore, na.rm=TRUE)
PMLclusdata_MockStdev<-sd(PMLclusdata$mock_avgzscore, na.rm=TRUE)
PMLclusdata_WtMean<-mean(PMLclusdata$hsv_avgzscore, na.rm=TRUE)
PMLclusdata_WtStdev<-sd(PMLclusdata$hsv_avgzscore, na.rm=TRUE)
PMLclusdata_MutMean<-mean(PMLclusdata$dICP0_avgzscore, na.rm=TRUE)
PMLclusdata_MutStdev<-sd(PMLclusdata$dICP0_avgzscore, na.rm=TRUE)
ClusData["PMLclus","MockMeanZscore"]<-PMLclusdata_MockMean
ClusData["PMLclus","MockStdevZscore"]<-PMLclusdata_MockStdev
ClusData["PMLclus","WtMeanZscore"]<-PMLclusdata_WtMean
ClusData["PMLclus","WtStdevZscore"]<-PMLclusdata_WtStdev
ClusData["PMLclus","MutMeanZscore"]<-PMLclusdata_MutMean
ClusData["PMLclus","MutStdevZscore"]<-PMLclusdata_MutStdev
PMLclusFCdata<-repfc[row.names(repfc) %in% PMLclus, ]
PMLclusFCdata_MockWtFcMean<-mean(PMLclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
PMLclusFCdata_MockWtFcStdev<-sd(PMLclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
PMLclusFCdata_MockMutFcMean<-mean(PMLclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
PMLclusFCdata_MockMutFcStdev<-sd(PMLclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
PMLclusFCdata_WtMutFcMean<-mean(PMLclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
PMLclusFCdata_WtMutFcStdev<-sd(PMLclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ClusData["PMLclus","MockWtMeanFC"]<-PMLclusFCdata_MockWtFcMean
ClusData["PMLclus","MockWtStdevFC"]<-PMLclusFCdata_MockWtFcStdev
ClusData["PMLclus","MockMutMeanFC"]<-PMLclusFCdata_MockMutFcMean
ClusData["PMLclus","MockMutStdevFC"]<-PMLclusFCdata_MockMutFcStdev
ClusData["PMLclus","WtMutMeanFC"]<-PMLclusFCdata_WtMutFcMean
ClusData["PMLclus","WtMutStdevFC"]<-PMLclusFCdata_WtMutFcStdev
#IFI16
IFI16clus<-vector()
IFI16clus<-scan("IFI16clus.txt", what='complex')
IFI16clusdata<-imputeddata_zscores[row.names(imputeddata_zscores) %in% IFI16clus, ]
IFI16clusdata_MockMean<-mean(IFI16clusdata$mock_avgzscore, na.rm=TRUE)
IFI16clusdata_MockStdev<-sd(IFI16clusdata$mock_avgzscore, na.rm=TRUE)
IFI16clusdata_WtMean<-mean(IFI16clusdata$hsv_avgzscore, na.rm=TRUE)
IFI16clusdata_WtStdev<-sd(IFI16clusdata$hsv_avgzscore, na.rm=TRUE)
IFI16clusdata_MutMean<-mean(IFI16clusdata$dICP0_avgzscore, na.rm=TRUE)
IFI16clusdata_MutStdev<-sd(IFI16clusdata$dICP0_avgzscore, na.rm=TRUE)
ClusData["IFI16clus","MockMeanZscore"]<-IFI16clusdata_MockMean
ClusData["IFI16clus","MockStdevZscore"]<-IFI16clusdata_MockStdev
ClusData["IFI16clus","WtMeanZscore"]<-IFI16clusdata_WtMean
ClusData["IFI16clus","WtStdevZscore"]<-IFI16clusdata_WtStdev
ClusData["IFI16clus","MutMeanZscore"]<-IFI16clusdata_MutMean
ClusData["IFI16clus","MutStdevZscore"]<-IFI16clusdata_MutStdev
IFI16clusFCdata<-repfc[row.names(repfc) %in% IFI16clus, ]
IFI16clusFCdata_MockWtFcMean<-mean(IFI16clusFCdata$mock_hsv_avgfc, na.rm=TRUE)
IFI16clusFCdata_MockWtFcStdev<-sd(IFI16clusFCdata$mock_hsv_avgfc, na.rm=TRUE)
IFI16clusFCdata_MockMutFcMean<-mean(IFI16clusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
IFI16clusFCdata_MockMutFcStdev<-sd(IFI16clusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
IFI16clusFCdata_WtMutFcMean<-mean(IFI16clusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
IFI16clusFCdata_WtMutFcStdev<-sd(IFI16clusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ClusData["IFI16clus","MockWtMeanFC"]<-IFI16clusFCdata_MockWtFcMean
ClusData["IFI16clus","MockWtStdevFC"]<-IFI16clusFCdata_MockWtFcStdev
ClusData["IFI16clus","MockMutMeanFC"]<-IFI16clusFCdata_MockMutFcMean
ClusData["IFI16clus","MockMutStdevFC"]<-IFI16clusFCdata_MockMutFcStdev
ClusData["IFI16clus","WtMutMeanFC"]<-IFI16clusFCdata_WtMutFcMean
ClusData["IFI16clus","WtMutStdevFC"]<-IFI16clusFCdata_WtMutFcStdev
#DNAPK
DNAPKclus<-vector()
DNAPKclus<-scan("DNAPKclus.txt", what='complex')
DNAPKclusdata<-imputeddata_zscores[row.names(imputeddata_zscores) %in% DNAPKclus, ]
DNAPKclusdata_MockMean<-mean(DNAPKclusdata$mock_avgzscore, na.rm=TRUE)
DNAPKclusdata_MockStdev<-sd(DNAPKclusdata$mock_avgzscore, na.rm=TRUE)
DNAPKclusdata_WtMean<-mean(DNAPKclusdata$hsv_avgzscore, na.rm=TRUE)
DNAPKclusdata_WtStdev<-sd(DNAPKclusdata$hsv_avgzscore, na.rm=TRUE)
DNAPKclusdata_MutMean<-mean(DNAPKclusdata$dICP0_avgzscore, na.rm=TRUE)
DNAPKclusdata_MutStdev<-sd(DNAPKclusdata$dICP0_avgzscore, na.rm=TRUE)
ClusData["DNAPKclus","MockMeanZscore"]<-DNAPKclusdata_MockMean
ClusData["DNAPKclus","MockStdevZscore"]<-DNAPKclusdata_MockStdev
ClusData["DNAPKclus","WtMeanZscore"]<-DNAPKclusdata_WtMean
ClusData["DNAPKclus","WtStdevZscore"]<-DNAPKclusdata_WtStdev
ClusData["DNAPKclus","MutMeanZscore"]<-DNAPKclusdata_MutMean
ClusData["DNAPKclus","MutStdevZscore"]<-DNAPKclusdata_MutStdev
DNAPKclusFCdata<-repfc[row.names(repfc) %in% DNAPKclus, ]
DNAPKclusFCdata_MockWtFcMean<-mean(DNAPKclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
DNAPKclusFCdata_MockWtFcStdev<-sd(DNAPKclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
DNAPKclusFCdata_MockMutFcMean<-mean(DNAPKclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
DNAPKclusFCdata_MockMutFcStdev<-sd(DNAPKclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
DNAPKclusFCdata_WtMutFcMean<-mean(DNAPKclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
DNAPKclusFCdata_WtMutFcStdev<-sd(DNAPKclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ClusData["DNAPKclus","MockWtMeanFC"]<-DNAPKclusFCdata_MockWtFcMean
ClusData["DNAPKclus","MockWtStdevFC"]<-DNAPKclusFCdata_MockWtFcStdev
ClusData["DNAPKclus","MockMutMeanFC"]<-DNAPKclusFCdata_MockMutFcMean
ClusData["DNAPKclus","MockMutStdevFC"]<-DNAPKclusFCdata_MockMutFcStdev
ClusData["DNAPKclus","WtMutMeanFC"]<-DNAPKclusFCdata_WtMutFcMean
ClusData["DNAPKclus","WtMutStdevFC"]<-DNAPKclusFCdata_WtMutFcStdev
#ATRX
ATRXclus<-vector()
ATRXclus<-scan("ATRXclus.txt", what='complex')
ATRXclusdata<-imputeddata_zscores[row.names(imputeddata_zscores) %in% ATRXclus, ]
ATRXclusdata_MockMean<-mean(ATRXclusdata$mock_avgzscore, na.rm=TRUE)
ATRXclusdata_MockStdev<-sd(ATRXclusdata$mock_avgzscore, na.rm=TRUE)
ATRXclusdata_WtMean<-mean(ATRXclusdata$hsv_avgzscore, na.rm=TRUE)
ATRXclusdata_WtStdev<-sd(ATRXclusdata$hsv_avgzscore, na.rm=TRUE)
ATRXclusdata_MutMean<-mean(ATRXclusdata$dICP0_avgzscore, na.rm=TRUE)
ATRXclusdata_MutStdev<-sd(ATRXclusdata$dICP0_avgzscore, na.rm=TRUE)
ClusData["ATRXclus","MockMeanZscore"]<-ATRXclusdata_MockMean
ClusData["ATRXclus","MockStdevZscore"]<-ATRXclusdata_MockStdev
ClusData["ATRXclus","WtMeanZscore"]<-ATRXclusdata_WtMean
ClusData["ATRXclus","WtStdevZscore"]<-ATRXclusdata_WtStdev
ClusData["ATRXclus","MutMeanZscore"]<-ATRXclusdata_MutMean
ClusData["ATRXclus","MutStdevZscore"]<-ATRXclusdata_MutStdev
ATRXclusFCdata<-repfc[row.names(repfc) %in% ATRXclus, ]
ATRXclusFCdata_MockWtFcMean<-mean(ATRXclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
ATRXclusFCdata_MockWtFcStdev<-sd(ATRXclusFCdata$mock_hsv_avgfc, na.rm=TRUE)
ATRXclusFCdata_MockMutFcMean<-mean(ATRXclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
ATRXclusFCdata_MockMutFcStdev<-sd(ATRXclusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
ATRXclusFCdata_WtMutFcMean<-mean(ATRXclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ATRXclusFCdata_WtMutFcStdev<-sd(ATRXclusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ClusData["ATRXclus","MockWtMeanFC"]<-ATRXclusFCdata_MockWtFcMean
ClusData["ATRXclus","MockWtStdevFC"]<-ATRXclusFCdata_MockWtFcStdev
ClusData["ATRXclus","MockMutMeanFC"]<-ATRXclusFCdata_MockMutFcMean
ClusData["ATRXclus","MockMutStdevFC"]<-ATRXclusFCdata_MockMutFcStdev
ClusData["ATRXclus","WtMutMeanFC"]<-ATRXclusFCdata_WtMutFcMean
ClusData["ATRXclus","WtMutStdevFC"]<-ATRXclusFCdata_WtMutFcStdev
#SLFN5
SLFN5clus<-vector()
SLFN5clus<-scan("SLFN5clus.txt", what='complex')
SLFN5clusdata<-imputeddata_zscores[row.names(imputeddata_zscores) %in% SLFN5clus, ]
SLFN5clusdata_MockMean<-mean(SLFN5clusdata$mock_avgzscore, na.rm=TRUE)
SLFN5clusdata_MockStdev<-sd(SLFN5clusdata$mock_avgzscore, na.rm=TRUE)
SLFN5clusdata_WtMean<-mean(SLFN5clusdata$hsv_avgzscore, na.rm=TRUE)
SLFN5clusdata_WtStdev<-sd(SLFN5clusdata$hsv_avgzscore, na.rm=TRUE)
SLFN5clusdata_MutMean<-mean(SLFN5clusdata$dICP0_avgzscore, na.rm=TRUE)
SLFN5clusdata_MutStdev<-sd(SLFN5clusdata$dICP0_avgzscore, na.rm=TRUE)
ClusData["SLFN5clus","MockMeanZscore"]<-SLFN5clusdata_MockMean
ClusData["SLFN5clus","MockStdevZscore"]<-SLFN5clusdata_MockStdev
ClusData["SLFN5clus","WtMeanZscore"]<-SLFN5clusdata_WtMean
ClusData["SLFN5clus","WtStdevZscore"]<-SLFN5clusdata_WtStdev
ClusData["SLFN5clus","MutMeanZscore"]<-SLFN5clusdata_MutMean
ClusData["SLFN5clus","MutStdevZscore"]<-SLFN5clusdata_MutStdev
SLFN5clusFCdata<-repfc[row.names(repfc) %in% SLFN5clus, ]
SLFN5clusFCdata_MockWtFcMean<-mean(SLFN5clusFCdata$mock_hsv_avgfc, na.rm=TRUE)
SLFN5clusFCdata_MockWtFcStdev<-sd(SLFN5clusFCdata$mock_hsv_avgfc, na.rm=TRUE)
SLFN5clusFCdata_MockMutFcMean<-mean(SLFN5clusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
SLFN5clusFCdata_MockMutFcStdev<-sd(SLFN5clusFCdata$mock_dICP0_avgfc, na.rm=TRUE)
SLFN5clusFCdata_WtMutFcMean<-mean(SLFN5clusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
SLFN5clusFCdata_WtMutFcStdev<-sd(SLFN5clusFCdata$hsv_dICP0_avgfc, na.rm=TRUE)
ClusData["SLFN5clus","MockWtMeanFC"]<-SLFN5clusFCdata_MockWtFcMean
ClusData["SLFN5clus","MockWtStdevFC"]<-SLFN5clusFCdata_MockWtFcStdev
ClusData["SLFN5clus","MockMutMeanFC"]<-SLFN5clusFCdata_MockMutFcMean
ClusData["SLFN5clus","MockMutStdevFC"]<-SLFN5clusFCdata_MockMutFcStdev
ClusData["SLFN5clus","WtMutMeanFC"]<-SLFN5clusFCdata_WtMutFcMean
ClusData["SLFN5clus","WtMutStdevFC"]<-SLFN5clusFCdata_WtMutFcStdev

###write.table(ClusData, file="PCACLusterDataAnalysis_ZscoreFoldChange.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)

}


