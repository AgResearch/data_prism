# 
#------------------------------------------------------------------------- 
# References : 
# http://www.inside-r.org/packages/cran/gplots/docs/heatmap.2
# http://stackoverflow.com/questions/21427863/how-to-add-more-margin-to-a-heatmap-2-plot-with-the-png-device
#-------------------------------------------------------------------------

# libraries
library(RColorBrewer)
setwd("/dataset/bioinformatics_dev/active/qb_paper/qb_paper/figure_1")

load_datamatrix2 <- function() {
   datamatrix2<<-read.table("datamatrix2.dat", header=TRUE, row.names=1, sep="\t")
}

build_datamatrix2 <- function() {
   # this resamples the original data

   datamatrix<-read.table("est_contig_refseq.txt", header=TRUE, row.names=1, sep="\t")
   datamatrix1<-datamatrix[,2:13]


   # initialise a matrix of selected genes - note , a random sample 
   datamatrix2<-datamatrix1[grep("nothing",datamatrix1$gene_name, ignore.case = TRUE),]
   datamatrix2<-datamatrix2[sample(nrow(datamatrix2), 0),]

   genes<-datamatrix1[grep("^PRL",datamatrix1$gene_name, ignore.case = TRUE),]
   genes<-genes[sample(nrow(genes), 20),]
   datamatrix2<-rbind(datamatrix2,genes)
   genes<-datamatrix1[grep("Sep15",datamatrix1$gene_name, ignore.case = TRUE),]
   genes<-genes[sample(nrow(genes), 20),]
   datamatrix2<-rbind(datamatrix2,genes)
   genes<-datamatrix1[grep("EEF1A1",datamatrix1$gene_name, ignore.case = TRUE),]
   genes<-genes[sample(nrow(genes), 20),]
   datamatrix2<-rbind(datamatrix2,genes)
   genes<-datamatrix1[grep("keratin 4",datamatrix1$gene_name, ignore.case = TRUE),]
   genes<-genes[sample(nrow(genes), 20),]
   datamatrix2<-rbind(datamatrix2,genes)
   genes<-datamatrix1[grep("keratin 8",datamatrix1$gene_name, ignore.case = TRUE),]
   genes<-genes[sample(nrow(genes), 20),]
   datamatrix2<-rbind(datamatrix2,genes)
   genes<-datamatrix1[grep("csn1s1",datamatrix1$gene_name, ignore.case = TRUE),]
   genes<-genes[sample(nrow(genes), 20),]
   datamatrix2<<-rbind(datamatrix2,genes)
}

save_datamatrix2 <- function() {
   # save datamatrix2 to a table which we can load 
   write.table(datamatrix2, file=paste("datamatrix2.dat"),row.names=TRUE,sep="\t")
}


do_plot <- function() {
   # for each gene in the gene_config table, use the regexp 
   # in that table to locate the matching ESTs, and hence set the 
   # appropriate plot colour 
   gene_config_table<-read.table("gene_config.txt", header=TRUE, row.names=1, sep="\t")
   point_colours <- rep(NA, nrow(datamatrix2))
   point_symbols <- rep('?', nrow(datamatrix2))
   point_count <- 1
   for(rowname in rownames(gene_config_table)) {
      regexp <-  toString(gene_config_table[rowname, "regexp"])
      colour <- paste("#", toString(gene_config_table[rowname, "colour"]),sep="")
      symbol <- toString(gene_config_table[rowname, "symbol"]) 
      point_colours[ grep(regexp,datamatrix2$gene, ignore.case = TRUE) ] <- colour
      point_symbols[ grep(regexp,datamatrix2$gene, ignore.case = TRUE) ] <- symbol
   }

   # euclidean 
   distances = dist(datamatrix2[,1:10])

   p1<- princomp(datamatrix2[,1:10])
   jpeg(filename = "figure_1.jpg", width=800, height=800) 
   plot.default(p1$scores[,0:2],col=point_colours, pch=point_symbols, cex=2.0,
       xlab="First Principal Component", ylab="Second Principal Component", cex.axis=1.5, cex.lab=1.5)

   title("PCA based metric structure for selected sequence spectra", cex.main=2.0)
   dev.off()
}

load_datamatrix2()
#build_datamatrix2()
do_plot()
#save_datamatrix2()




