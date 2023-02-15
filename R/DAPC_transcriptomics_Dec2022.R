#' Discriminant analysis of principal components (DAPC) of transcriptomics dataset
#'
#' @import adegenet
#'
#' @return
#' @export

DAPC_transcriptomics_Dec2022 = function(){
  fpath = system.file("extdata", "SF2_Circadian_PX_transcriptome.csv",
                      package = "PublicationPleuromammaOmics")

  #read in transcriptomics data matrix
  trans.dat<-read.csv(fpath, header=T, row.names=1)
  #data starts in column 6

  #sum rows containing data
  trans.dat$sum<-rowSums(trans.dat)

  #keep only rows that sum to non-zero values
  trans.dat2<-subset(trans.dat, sum>0, select=TS5b:TS8L)

  #check data for missing values and outliers
  gsg = WGCNA::goodSamplesGenes(trans.dat2, verbose=3);
  gsg$allOK
  #response is TRUE so no genes or samples need to be removed

  #this custom theme will be used for the DAPC plot
  theme_custom <- function() {
    theme_bw(base_size = 10) %+replace%    #, base_family = "Arial"
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.background = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
      )
  }


  #timepoint metadata following order of column names
  trans.grp<-c(rep("22:00",6), rep('02:00',6), rep('09:00', 6), rep('15:00',5))

  #determine how many PCs to keep for the analysis
  dapc.test<-dapc(t(trans.dat2), trans.grp, n.da=100, n.pca=50)
  ascore.test<-optim.a.score(dapc.test, n.sim=5)

  #optimal number of PCs is 11
  xval<-xvalDapc(t(trans.dat2), trans.grp, n.pca.max=300, training.set=0.9,
                 result='groupMean', center=T, scale=F, n.pca=NULL, n.rep=30,
                 xval.plot=T)

  #DAPC analysis
  dapc.trans<-dapc(t(trans.dat2), trans.grp, n.da=100, n.pca=11)

  dapc<-tibble(sample=rownames(dapc.trans$ind.coord),
               grp=dapc.trans$grp,
               LD1=dapc.trans$ind.coord[,1],
               LD2=dapc.trans$ind.coord[,2])
  dapc<-dapc %>%
    group_by(grp) %>%
    summarize(c1=mean(LD1),
              c2=mean(LD2)) %>%
    full_join(dapc)

  varexpl<-round((dapc.trans$eig/sum(dapc.trans$eig))[1:2] * 100, 1)

  if(!file.exists("figs"))
    dir.create("figs")

  #create plot of DAPC
  dapc.fig<-
    ggplot(dapc, aes(shape = factor(grp))) +
    geom_segment(mapping=aes(x=LD1, y=LD2, xend=c1, yend=c2), lwd=0.25, col='grey') +
    geom_point(aes(x=c1, y=c2, colour=factor(grp), fill=factor(grp)), size=4) +
    geom_point(aes(x=LD1, y=LD2, colour=factor(grp), fill=factor(grp)), size=1) +
    labs(x=paste0('LD1 [', varexpl[1], '%]'), y=paste0('LD2 [', varexpl[2], '%]')) +
    scale_shape_manual(name='timepoint',
                       labels=c('02:00', '09:00', '15:00', '22:00'),
                       values=c(rep(21,4))) +
    scale_colour_manual(name='timepoint',
                        labels=c('02:00', '09:00', '15:00', '22:00'),
                        values=c(rep('black',4))) +
    scale_fill_manual(name='timepoint',
                      labels=c('02:00', '09:00', '15:00', '22:00'),
                      values=c('blue', 'goldenrod1', 'yellow', "slategray3")) +
    theme_custom() +
    theme(plot.margin = margin(t=0.5, l=0.5, r=0.5, b=1, unit='cm')) +
    ggtitle('Transcriptomics') +
    theme(plot.title=element_text(hjust = 0.5, size=16))

  pdf("figs/DAPC transcriptomics Dec2022.pdf", width=7, height=7, pointsize=8)

  print(dapc.fig)

  dev.off()

  #explort the LD loadings for each axis
  set.seed(983)
  a<-data.frame(dapc.trans$var.contr)
  gene_rank<-c(1:nrow(a))

  if(!file.exists("out"))
    dir.create("out")

  #LD1
  ldl1<-data.frame(a[c(1)])
  ld1<-data.frame(ldl1[order(-ldl1$LD1),,drop=F])
  gene.ld1<-ld1[c(1:60),,drop=F]
  write.table(gene.ld1, 'out/DAPC transcriptomics LD1 proteins Dec2022.txt', quote=F)

  #LD2
  ldl2 <-data.frame(a[c(2)])
  ld2<-data.frame(ldl2[order(-ldl2$LD2),,drop=F])
  gene.ld2<-ld2[c(1:60),,drop=F]
  write.table(gene.ld2, 'out/DAPC transcriptomics LD2 proteins Dec2022.txt', quote=F)
}
