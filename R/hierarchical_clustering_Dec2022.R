#' Clustering
#'
#' @import vegan
#' @import ggthemes
#' @import reshape
#' @import ggplot2
#'
#' @return
#' @export

hierarchical_clustering_Dec2022 = function(){
  fpath = system.file("extdata", "SF2_Circadian_PX_transcriptome.csv",
                      package = "PublicationPleuromammaOmics")
  trans.dat<-read.csv(fpath, header=T, row.names=1)

  #keep only rows that sum to non-zero values
  trans.dat2<-subset(trans.dat, select=TS5b:TS8L)
  trans.dat2$sum<-rowSums(trans.dat2)
  trans.dat2<-subset(trans.dat2, sum>0, select=TS5b:TS8L)

  #average by time point
  TS22<-cbind(trans.dat2$TS5b, trans.dat2$TS5c, trans.dat2$TS5d,
              trans.dat2$TS9b, trans.dat2$TS9d, trans.dat2$TS8b)
  TS2<-cbind(trans.dat2$TS6a, trans.dat2$TS6b, trans.dat2$TS6c,
             trans.dat2$TS6d, trans.dat2$TS6e,trans.dat2$TS6f)
  TS9<-cbind(trans.dat2$TS7a, trans.dat2$TS7c, trans.dat2$TS7e,
             trans.dat2$TS7f, trans.dat2$TS7g, trans.dat2$TS7h)
  TS15<-cbind(trans.dat2$TS8a, trans.dat2$TS8d, trans.dat2$TS8J,
              trans.dat2$TS8K, trans.dat2$TS8L)

  TS22.avg<-rowMeans(TS22)
  TS2.avg<-rowMeans(TS2)
  TS9.avg<-rowMeans(TS9)
  TS15.avg<-rowMeans(TS15)

  trans.avg.df<-data.frame(TS22.avg, TS2.avg, TS9.avg, TS15.avg)
  rownames(trans.avg.df)<-rownames(trans.dat2)

  #create a bray-curtis dissimilarity matrix
  trans.bray<-vegdist(trans.avg.df, method='bray')

  #cluster matrix using the average clustering method
  trans.clust<-hclust(trans.bray, method='average')

  #plot cluster dendrogram without labels
  plot(trans.clust, labels=F)

  #plot rectangles with height cutoff of 0.4
  rect.hclust(trans.clust, h=0.4)

  #cut the dendrogram at h=0.4
  trans0.4cut<-cutree(trans.clust, h=0.4)

  if(!file.exists("out"))
    dir.create("out")

  #write cluster assignments to a csv file
  write.csv(trans0.4cut,
            'out/cluster assignments transcriptomics h 0.4 Dec2022.csv', quote=F)

  #this is clunky: in excel, change cluster assignment csv file headings to
  #"Gene" and "cluster"

  #read the file back in with new column names
  clust0.4<-read.csv('out/cluster assignments transcriptomics h 0.4 Dec2022.csv',
                     header=T, row.names=1)
  colnames(clust0.4) = c("cluster")

  #merge cluster assignments with protein abundance values
  clust.genexp<-merge(x=clust0.4, y=trans.avg.df, by='row.names', all.x=T)

  #melt the data set to reformat for plotting
  melt.trans0.4<-melt(clust.genexp, id.vars=c('Row.names', 'cluster'))

  #make the plot: each cluster will be its own plot in the grid of plots.
  #Each plot is protein abundance plotted against time point where each line in
  #the plot represents a single protein's abundance. Also uses geom_smooth to
  #add loess smoothed curve to the plots

  if(!file.exists("figs"))
    dir.create("figs")

  g<-
    ggplot(melt.trans0.4, aes(x=variable, y=value, group=Row.names)) +
    geom_line(alpha=0.1) +
    theme_bw() +
    facet_wrap(~cluster, scales='free_y') +
    labs(x='Time Point', y='Averaged Gene Expression') +
    theme(axis.text.x = element_text(angle=90,size=6)) +
    geom_smooth(aes(group=1),method="loess", se=FALSE, span=0.6) +
    annotate("rect", xmin=0, xmax=2.5, ymin=0, ymax=Inf, alpha=0.3, fill='blue') +
    annotate("rect", xmin=2.5, xmax=4.5, ymin=0, ymax=Inf, alpha=0.3, fill='yellow') +
    scale_x_discrete('Time point',labels=c('22:00', '02:00', '09:00', '15:00'))

  png("figs/transcriptomics clustering Dec2022.png", width=15, height=15,
      units = "in", res = 150)

  print(g)

  dev.off()
}
