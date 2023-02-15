#' Proteomics NMDS and hierarchical clustering
#'
#' @import vegan
#'
#' @return
#' @export

AE1910_metazoan_Dec2022 = function(){

  fpath = system.file("extdata", "ABACUS_copepods_metazoan_output.csv",
                      package = "PublicationPleuromammaOmics")

  cop.dat<-read.csv(fpath, header=T)

  #keep only entries that are predicted copepod proteins
  #make new column in cop.dat that includes protein name before "|" - 24 characters
  cop.dat$TranscriptID <- substr(cop.dat$PROTID, 1, 24)
  #merge with list of protein names from metazoan transcriptome fasta and keep only those that match
  fpath = system.file("extdata", "Metazoan_isoforms.txt",
                      package = "PublicationPleuromammaOmics")
  trans.list<-read.table(fpath, header=F, sep='\t')
  names(trans.list)[names(trans.list)=='V1']<-'TranscriptID'

  only.cop<-merge(x=cop.dat, y=trans.list, by='TranscriptID')
  #2551 proteins (original file had 2767)

  #add PROTID column as row names
  row.names(only.cop)<-only.cop$PROTID

  adjnsaf<-select(only.cop, contains('ADJNSAF'))

  #Keep only proteins with at least 2 unique peptides
  nsaf.uniq<-cbind(adjnsaf, only.cop$ALL_NUMPEPSUNIQ)
  twopeps<-subset(nsaf.uniq,
                  select= X20190815_UWPRQE_COPEPODS_01_ADJNSAF:X20190815_UWPRQE_COPEPODS_77_ADJNSAF,
                  nsaf.uniq[,77]>1)

  #NMDS
  #Remove irrelevant samples
  cop.sub1<-subset(twopeps, select=01:61)
  cop.TS<-subset(cop.sub1, select=-13)

  cop.t<-t(cop.TS)
  cop.tra<-(cop.t+1)
  cop.tra<-data.trans(cop.tra, method='log', plot=F)

  nmds.techreps<-metaMDS(cop.tra, distance='bray', k=2, trymax=100,
                         autotransform=F)

  ordiplot(nmds.techreps, choices=c(1,2), type='text', display='sites', cex=0.5)

  #remove samples 2-4 that do not cluster well with replicates
  cop.TS2<-subset(cop.TS, select=c(01,05:60))

  #average technical replicates
  TS5.1<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_16_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_22_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_58_ADJNSAF)
  TS5.2<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_08_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_32_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_57_ADJNSAF)
  TS5.3<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_01_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_27_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_56_ADJNSAF)
  TS5.4<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_05_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_23_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_46_ADJNSAF)
  TS6.1<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_13_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_26_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_54_ADJNSAF)
  TS6.2<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_28_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_35_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_51_ADJNSAF)
  TS6.3<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_03_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_21_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_43_ADJNSAF)
  TS6.4<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_07_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_33_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_50_ADJNSAF)
  TS7.1<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_29_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_36_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_49_ADJNSAF)
  TS7.2<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_17_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_40_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_42_ADJNSAF)
  TS7.3<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_12_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_31_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_59_ADJNSAF)
  TS7.4<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_14_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_38_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_52_ADJNSAF)
  TS8.1<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_11_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_20_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_60_ADJNSAF)
  TS8.2<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_02_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_34_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_48_ADJNSAF)
  TS8.3<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_30_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_41_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_45_ADJNSAF)
  TS8.4<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_04_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_39_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_44_ADJNSAF)
  TS9.1<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_06_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_37_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_61_ADJNSAF)
  TS9.2<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_09_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_19_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_55_ADJNSAF)
  TS9.3<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_10_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_25_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_53_ADJNSAF)
  TS9.4<-cbind(cop.TS2$X20190815_UWPRQE_COPEPODS_15_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_24_ADJNSAF,
               cop.TS2$X20190815_UWPRQE_COPEPODS_47_ADJNSAF)

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
  rownames(TS.df)<-rownames(cop.TS2)

  TS.t<-t(TS.df)
  TS.tra<-(TS.t+1)
  TS.tra<-data.trans(TS.tra, method='log', plot=F)

  nmds.TS<-metaMDS(TS.tra, distance='bray', k=2, trymax=100, autotransform=F)

  TS.eigen<-envfit(nmds.TS$points, TS.tra, perm=1000)

  if(!file.exists("figs"))
    dir.create("figs")

  pdf('figs/Pxiphias proteomics NMDS Jan2023.pdf', width=7, height=7,
      pointsize=8)
  ordiplot(nmds.TS, choices=c(1,2), type='text', display='sites', cex=0.5)
  TS.grp<-c(rep('22:00',4), rep('02:00',4), rep('09:00',4), rep('15:00',4),
            rep('22:00',4))
  fig.TS1<-ordiplot(nmds.TS, choices=c(1,2), type='none', display='sites')
  points(fig.TS1, 'sites', bg=c(rep("slategray3",4), rep("blue",4),
                                rep("goldenrod1",4), rep("yellow",4),
                                rep("slategray3",4)), pch=21, cex=1.2)
  ordihull(fig.TS1,groups=TS.grp,draw='lines',col='grey75',label=F)
  legend('topleft', legend=c('22:00', '02:00', '09:00', '15:00'), pch=19,
         col=c('slategray3', 'blue', 'goldenrod1', 'yellow'))
  dev.off()

  #ANOSIM day vs night
  dayvsnight<-c(rep('night',8), rep('day',8), rep('night',4))

  TS.row<-data.stand(TS.t, method='total', margin='row', plot=F)
  TS.d<-vegdist(TS.row, 'bray')
  dayvsnight.anosim<-anosim(TS.d, grouping=dayvsnight)
  summary(dayvsnight.anosim)


  #ANOSIM time points
  time.points.anosim<-anosim(TS.d, grouping=TS.grp)
  summary(time.points.anosim)

  #Clustering
  TS5<-cbind(TS.df$TS5.1.avg, TS.df$TS5.2.avg, TS.df$TS5.3.avg, TS.df$TS5.4.avg)
  TS6<-cbind(TS.df$TS6.1.avg, TS.df$TS6.2.avg, TS.df$TS6.3.avg, TS.df$TS6.4.avg)
  TS7<-cbind(TS.df$TS7.1.avg, TS.df$TS7.2.avg, TS.df$TS7.3.avg, TS.df$TS7.4.avg)
  TS8<-cbind(TS.df$TS8.1.avg, TS.df$TS8.2.avg, TS.df$TS8.3.avg, TS.df$TS8.4.avg)
  TS9<-cbind(TS.df$TS9.1.avg, TS.df$TS9.2.avg, TS.df$TS9.3.avg, TS.df$TS9.4.avg)

  TS22.avg<-rowMeans(TS5, TS9)
  TS2.avg<-rowMeans(TS6)
  TS9.avg<-rowMeans(TS7)
  TS15.avg<-rowMeans(TS8)

  TS.avg.df<-data.frame(TS22.avg, TS2.avg, TS9.avg, TS15.avg)
  rownames(TS.avg.df)<-rownames(cop.TS2)

  TS.avg.df$sum<-rowSums(TS.avg.df)
  TS.avg.df2<-subset(TS.avg.df, sum>0, select=-sum)
  TS.avg.bray<-vegdist(TS.avg.df2, method='bray')
  TS.avg.clust<-hclust(TS.avg.bray, method='average')
  plot(TS.avg.clust, labels=F)
  TS.cut<-cutree(TS.avg.clust, h=0.5)

  if(!file.exists("out"))
    dir.create("out")

  write.csv(TS.cut, 'out/cluster assignments TS Jan2023.csv', quote=F)

  prot.clust.TS2<-read.csv('out/cluster assignments TS Jan2023.csv', header=T,
                           row.names=1)
  clust.nsaf.TS2<-merge(x=prot.clust.TS2, y=TS.avg.df2, by='row.names', all.x=T)

  library(reshape2)
  melt.TS2<-melt(clust.nsaf.TS2, id.vars=c('Row.names', 'x'))

  library(ggplot2)
  gg.TS<-
    ggplot(melt.TS2, aes(x=variable, y=value, group=Row.names)) +
    geom_line(alpha=0.1) +
    theme_bw() +
    facet_wrap(~x, scales='free_y') +
    labs(x='Time Point', y='Averaged Normalized Spectral Abundance Factor') +
    theme(axis.text.x = element_text(angle=90,size=6)) +
    geom_smooth(aes(group=1),method="loess", se=FALSE, span=0.6) +
    annotate("rect", xmin=0, xmax=2.5, ymin=0, ymax=Inf, alpha=0.3,
             fill='blue') +
    annotate("rect", xmin=2.5, xmax=4.5, ymin=0, ymax=Inf, alpha=0.3,
             fill='yellow') +
    scale_x_discrete('Time point',labels=c('22:00', '02:00', '09:00', '15:00'))

  png('figs/Pxiphias cluster plot Jan2023.png', width=7, height=7, units='in',
      res=150)
  print(gg.TS)
  dev.off()
}
