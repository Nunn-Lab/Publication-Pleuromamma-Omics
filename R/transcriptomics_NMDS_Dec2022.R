#' Nonmetric multidimensional scaling (NMDS) of transcriptomics dataset
#'
#' @return
#' @export

transcriptomics_NMDS_Dec2022 = function(){
  fpath = system.file("extdata", "SF2_Circadian_PX_transcriptome.csv",
                      package = "PublicationPleuromammaOmics")

  #read in transcriptomics data matrix
  trans.dat<-read.csv(fpath, header=T, row.names=1)

  #sum rows containing data
  trans.dat$sum<-rowSums(trans.dat)

  #keep only rows that sum to non-zero values
  trans.dat2<-subset(trans.dat, sum>0, select=TS5b:TS8L)
  #9a is an outlier

  #NMDS
  #transpose the data
  trans.t<-t(trans.dat2)

  #data transformation: Log(x+1)
  trans.tra<-data.trans((trans.t+1), method='log', plot=F)

  #do NMDS ordination
  trans.nmds<-metaMDS(trans.tra, distance='bray', k=2, trymax=100,
                      autotransform=F)

  #generate plot with column names
  ordiplot(trans.nmds, choices=c(1,2), type='text', display='sites', cex=0.5)

  #generate plot with colors & shapes to represent time points and time of day

  #Final NMDS figure
  trans.grp<-c(rep("22:00",6), rep('02:00',6), rep('09:00', 6), rep('15:00',5))

  if(!file.exists("figs"))
    dir.create("figs")

  pdf("figs/NMDS transcriptomics Dec2022.pdf", width=7, height=7, pointsize=8)
  fig.trans<-ordiplot(trans.nmds, choices=c(1,2), type='none', display='sites')
  points(fig.trans, 'sites', bg=c(rep('slategray3',6), rep('blue',6),
                                  rep('goldenrod1', 6), rep('yellow',5)),
         pch=21, cex=5)
  ordihull(fig.trans,groups=trans.grp,draw='lines',col='grey75',label=F)
  legend('topleft', legend=c('22:00', '02:00', '09:00', '15:00'), pch=19,
         col=c("slategray3", 'blue', 'goldenrod1', 'yellow'))
  dev.off()
}
