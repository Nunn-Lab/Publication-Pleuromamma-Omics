#' DAPC of copepod proteomics
#'
#' @import adegenet
#' @import tidyverse
#' @import WGCNA
#' @import ggplot2
#'
#' @return
#' @export

dapc_copepod_dda_Jan2023 = function(){

  #read in your data file
  fpath = system.file("extdata", "Copepod Metazoan Proteins NSAF.csv",
                      package = "PublicationPleuromammaOmics")
  copDDA<-read.csv(fpath, row.names=1)
  cop.TS<-subset(copDDA, select=01:61)
  cop.TS<-subset(cop.TS, select=-c(2,3,4,13))

  #Transpose the data frame and fix it up
  copDDA.t<-as.data.frame(t(cop.TS))

  #check data for missing values and outliers
  gsg = goodSamplesGenes(copDDA.t, verbose=3);
  gsg$allOK
  #if TRUE,go to line 34
  #if FALSE< some proteins have too many missing values and run the following code
  if (!gsg$allOK)
  {
    if(sum(!gsg$goodGenes)>0)
      dynamicTreeCut::printFlush(paste("Removing genes:",
                       paste(names(copDDA.t) [!gsg$goodGenes], collapse=", ")));
    if(sum(!gsg$goodSamples)>0)
      dynamicTreeCut::printFlush(paste("Removing samples:",
                       paste(rownames(copDDA.t) [!gsg$goodSamples], collapse=", ")));
    copDDA.t = copDDA.t[gsg$goodSamples, gsg$goodGenes]
  }

  #average tech reps
  copDDA2<-data.frame(t(copDDA.t))
  TS5.1<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_16_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_22_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_58_ADJNSAF)
  TS5.2<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_08_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_32_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_57_ADJNSAF)
  TS5.3<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_01_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_27_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_56_ADJNSAF)
  TS5.4<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_05_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_23_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_46_ADJNSAF)
  TS6.1<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_26_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_54_ADJNSAF)
  TS6.2<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_28_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_35_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_51_ADJNSAF)
  TS6.3<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_21_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_43_ADJNSAF)
  TS6.4<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_07_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_33_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_50_ADJNSAF)
  TS7.1<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_29_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_36_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_49_ADJNSAF)
  TS7.2<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_17_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_40_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_42_ADJNSAF)
  TS7.3<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_12_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_31_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_59_ADJNSAF)
  TS7.4<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_14_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_38_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_52_ADJNSAF)
  TS8.1<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_11_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_20_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_60_ADJNSAF)
  TS8.2<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_34_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_48_ADJNSAF)
  TS8.3<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_30_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_41_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_45_ADJNSAF)
  TS8.4<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_39_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_44_ADJNSAF)
  TS9.1<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_06_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_37_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_61_ADJNSAF)
  TS9.2<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_09_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_19_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_55_ADJNSAF)
  TS9.3<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_10_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_25_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_53_ADJNSAF)
  TS9.4<-cbind(copDDA2$X20190815_UWPRQE_COPEPODS_15_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_24_ADJNSAF,
               copDDA2$X20190815_UWPRQE_COPEPODS_47_ADJNSAF)

  #run line of code below for each of the combined time points above
  TS5.1.avg<-rowMeans(TS5.1)
  TS.avg = list(rowMeans(TS5.1),rowMeans(TS5.2),rowMeans(TS5.3),rowMeans(TS5.4),
                rowMeans(TS6.1),rowMeans(TS6.2),rowMeans(TS6.3),rowMeans(TS6.4),
                rowMeans(TS7.1),rowMeans(TS7.2),rowMeans(TS7.3),rowMeans(TS7.4),
                rowMeans(TS8.1),rowMeans(TS8.2),rowMeans(TS8.3),rowMeans(TS8.4),
                rowMeans(TS9.1),rowMeans(TS9.2),rowMeans(TS9.3),rowMeans(TS9.4))
  names(TS.avg) = c("TS5.1.avg", "TS5.2.avg", "TS5.3.avg", "TS5.4.avg",
                    "TS6.1.avg", "TS6.2.avg", "TS6.3.avg", "TS6.4.avg",
                    "TS7.1.avg", "TS7.2.avg", "TS7.3.avg", "TS7.4.avg",
                    "TS8.1.avg", "TS8.2.avg", "TS8.3.avg", "TS8.4.avg",
                    "TS9.1.avg", "TS9.2.avg", "TS9.3.avg", "TS9.4.avg")

  TS.df<-data.frame(TS.avg[[1]], TS.avg[[2]], TS.avg[[3]], TS.avg[[4]],
                    TS.avg[[5]], TS.avg[[6]], TS.avg[[7]], TS.avg[[8]],
                    TS.avg[[9]], TS.avg[[10]], TS.avg[[11]], TS.avg[[12]],
                    TS.avg[[13]], TS.avg[[14]], TS.avg[[15]], TS.avg[[16]],
                    TS.avg[[17]], TS.avg[[18]], TS.avg[[19]], TS.avg[[20]])
  colnames(TS.df) = names(TS.avg)
  rownames(TS.df)<-rownames(copDDA2)

  cop.DDAt2<-t(TS.df)

  theme_custom <- function() {
    theme_bw(base_size = 10) %+replace%    #, base_family = "Arial"
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.background = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(angle=45, hjust=1, vjust = 1)#,
        #legend.title = element_text(size = 8),
        #legend.text = element_text(size = 7)
      )
  }

  #timepoint metadata following order of colnames
  meta<-c('22:00', '22:00','22:00','22:00','02:00','02:00','02:00','02:00',
          '09:00','09:00','09:00','09:00', '15:00','15:00','15:00','15:00',
          '22:00','22:00','22:00','22:00')

  #how many PCs to keep
  dapc.test<-dapc(cop.DDAt2, meta, n.da=100, n.pca=50)
  ascore.test<-optim.a.score(dapc.test, n.sim=5)
  #keep 6 PCs

  xval<-xvalDapc(cop.DDAt2, meta, n.pca.max=300, training.set=0.9,
                 result='groupMean', center=T, scale=F, n.pca=NULL, n.rep=30,
                 xval.plot=T)
  #find x-axis value that corresponds to highest y-axis value;
  #choose x-axis number for n.pca below

  dapc.dda<-dapc(cop.DDAt2, meta, n.da=100, n.pca=6)
  ade4::scatter(dapc.dda, bg='white', scree.da=T, scree.pca=T, legend=T, solid=.4)

  dapc<-tibble(sample=rownames(dapc.dda$ind.coord),
               grp=dapc.dda$grp,
               LD1=dapc.dda$ind.coord[,1],
               LD2=dapc.dda$ind.coord[,2])
  dapc<-dapc %>%
    group_by(grp) %>%
    summarize(c1=mean(LD1),
              c2=mean(LD2)) %>%
    full_join(dapc)

  varexpl<-round((dapc.dda$eig/sum(dapc.dda$eig))[1:2] * 100, 1)

  #spider plot

  if(!file.exists("figs"))
    dir.create("figs")

  dapc.fig<-
    ggplot(dapc, aes(shape = factor(grp))) +
    geom_segment(mapping=aes(x=LD1, y=LD2, xend=c1, yend=c2),
                 lwd=0.25, col='grey') +
    geom_point(aes(x=c1, y=c2, fill=factor(grp)), size=4) +
    geom_point(aes(x=LD1, y=LD2, fill=factor(grp)), size=1) +
    labs(x=paste0('LD1 [', varexpl[1], '%]'),
         y=paste0('LD2 [', varexpl[2], '%]')) +
    scale_shape_manual(name='timepoint',
                       labels=c('02:00', '09:00', '15:00', '22:00'),
                       values=c(rep(21,4))) +
    scale_fill_manual(name='timepoint',
                      labels=c('02:00', '09:00', '15:00', '22:00'),
                      values=c( 'blue', 'goldenrod1', 'yellow', 'slategray3')) +
    theme_custom() +
    theme(plot.margin = margin(t=0.5, l=0.5, r=0.5, b=1, unit='cm')) +
    ggtitle('Proteomics') +
    theme(plot.title=element_text(hjust = 0.5, size=16))

  png('figs/Pxiphias DAPC proteomics Jan2023.png',
      width=7, height=7, units='in', res=150)
  print(dapc.fig)
  dev.off()


  #LD loadings
  set.seed(983)
  a<-data.frame(dapc.dda$var.contr)
  protein_rank<-c(1:nrow(a))

  if(!file.exists("out"))
    dir.create("out")

  #LD1
  ldl1<-data.frame(a[c(1)])
  ld1<-data.frame(ldl1[order(-ldl1$LD1),,drop=F])
  prot.ld1<-ld1[c(1:60),,drop=F]
  write.table(prot.ld1, 'out/DAPC DDA LD1 proteins Jan2023.txt', quote=F)
  #LD2
  ldl2 <-data.frame(a[c(2)])
  ld2<-data.frame(ldl2[order(-ldl2$LD2),,drop=F])
  prot.ld2<-ld2[c(1:60),,drop=F]
  write.table(prot.ld2, 'out/DAPC DDA LD2 proteins Jan2023.txt', quote=F)
}
