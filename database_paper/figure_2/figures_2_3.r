library(Heatplus)
library(RColorBrewer)
library("gplots")
setwd("/dataset/bioinformatics_dev/active/qb_paper/qb_paper/figure_2")


num_clust=300


build_clustered_data <- function() {

   datamatrix <<-read.table("annotated_data_matrix.txt", header=TRUE, sep="\t")  

 
   clustering <<- kmeans(datamatrix[,2:11], num_clust, iter.max=500)

   tcenters <<-t(clustering$centers)
   dists <<- dist(tcenters)

   h<-hclust(dists)
   dgram <<- as.dendrogram(h)

   mycolv<<-dgram
}

annotate_data <- function() {
#estname bgisheep_agbovine.csv   BovineVelevetSE_agbovine.csv    btau42_agbovine.csv     btau461_agbovine.csv    cs34_agbovine.csv       cs39_agbovine.csv       dfciBt_agbovine.csv dfciOa_agbovine.csv     umd2_agbovine.csv       umd3_agbovine.csv
#000110BPBA008196HT      11.541432341    18.4000872438   12.2328052503   11.3991443849   12.3656722052   13.5769074721   13.5761816869   13.2523143241   12.4262647547   12.1621290349
#000110BPBA008197HT      17.4956286514   18.4000872438   17.5182074691   17.4866072262   17.5751255708   18.1618699729   17.5761816869   18.2523143241   17.5555477716   17.5544464577

#intrepid$ head est_genes.txt
#estname accession       gene_name
#000110BPBA008196HT      XM_005223937.3  MBP
#000110BPBA008198HT      XM_005223937.3  MBP

   datamatrix <-read.table("space_zero.txt", header=TRUE, sep="\t")  
   est_gene<-read.table("est_genes.txt", header=TRUE, row.names=1, sep="\t")
   annotated_data_matrix <- cbind(datamatrix, est_gene[as.character(datamatrix$estname),"gene_name"])
   colnames(annotated_data_matrix) <- c(colnames(datamatrix),"gene_name")
   write.table(annotated_data_matrix, "annotated_data_matrix.txt", sep="\t")
}

build_plot_data <- function() {


   reductionFactor = 40
   slice = 0

   reductionSelector = sequence(nrow(datamatrix))
   reductionSelector <- subset(reductionSelector, reductionSelector %% reductionFactor == slice)
   datamatrix<<-datamatrix[reductionSelector,]

   orderedcols <- colnames(as.matrix(dists))
   datamatrix <<- datamatrix[,c("estname", orderedcols, "gene_name")]


   labelInterval=50  # 1=label every probe, 2=label every 2nd probe etc (for readability if many probes on plot)

   rowLabels_est <<- as.vector(datamatrix$estname)
   rowLabels_gene <<- as.vector(datamatrix[,"gene_name"])
   #rowLabels <<- rowLabels_gene 


   rownames(datamatrix)<<-datamatrix$estname
   datamatrix<<-datamatrix[,2:11] 


   # replace the sample filenames with better names 
   genome_names<-read.table("genome_names.dat", header=TRUE, row.names=1, sep="\t")
   colnames(datamatrix) <<- as.vector(genome_names[colnames(datamatrix),"genome_name"])   



   my_palette <<- brewer.pal(11,"Spectral") # a diverging palette

   blankSelector <- sequence(length(rowLabels_gene))
   blankSelector <- subset(blankSelector, blankSelector %% labelInterval != 1) # e.g. will get (2,3, 5,6, 8,9, ..)
                                                                            # so we will only label rows 1,4,7,10,13 etc)
   jpeg(filename =paste("hm1_internal", slice, ".jpg",sep=""), width=800, height=1200) # with dendrograms


   # run the heatmap, just to obtain the clustering index - not the final plot
   hm_internal<-heatmap.2(as.matrix(datamatrix),  scale = "none", dendrogram = "col",  
       Colv = mycolv, 
       trace = "none",
       #trace = "none", breaks = 5/11*seq(0,11), 
       #trace = "none", breaks =  -2 + 4/11*seq(0,11), 
       col = my_palette , key=TRUE, density.info="none", 
       #keysize=1.0, margin=c(9,9), cexRow=1.2, cexCol=.9, 
       keysize=1.0, margin=c(9,9), cexRow=1.5, cexCol=.9, 
       #lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(1.0, 2.0, 0 ), lhei=c(.5, 3) )
       #lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(0.1, .5, 0), lhei=c(.5, 3) )
       lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(0.1, .6, 0), lhei=c(.5, 3.5) )
   dev.off()

   # now edit the re-ordered vector of row labels, obtained from the heatmap object, so that only 
   # every nth label on the final plot has a non-empty string
   indexSelector <- hm_internal$rowInd[length(hm_internal$rowInd):1]    
   indexSelector <- indexSelector[blankSelector]
   rowLabels_gene[indexSelector] <<- rep('',length(indexSelector))
   rowLabels_est[indexSelector] <<- rep('',length(indexSelector))

}

do_figure2 <- function() {


   jpeg(filename = "figure_2.jpg", width=830, height=1400) # with dendrograms

   hm_gene_labels<<-heatmap.2(as.matrix(datamatrix),  scale = "none", 
    dendrogram = "col",  
    Colv = mycolv,
    #trace="none",
    trace = "none", breaks =  5 + 13/11*seq(0,11), 
    col = my_palette , key=TRUE, density.info="none", 
    #keysize=1.0, margin=c(17,25), cexRow=1.5, cexCol=1.6, 
    #keysize=1.0, margin=c(25,20), cexRow=1.5, cexCol=1.6, 
    keysize=1.0, margin=c(25,20), cexRow=1.4, cexCol=2.0, 

    lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, .6, 0 ), lhei=c(.2, 3) , labRow = rowLabels_gene)


   dev.off()

   jpeg(filename = "figure_2_estlables.jpg", width=830, height=1400) # with dendrograms

   hm_est_labels<<-heatmap.2(as.matrix(datamatrix),  scale = "none", 
    dendrogram = "col",  
    Colv = mycolv,
    #trace="none",
    trace = "none", breaks =  5 + 13/11*seq(0,11), 
    col = my_palette , key=TRUE, density.info="none", 
    #keysize=1.0, margin=c(17,25), cexRow=1.5, cexCol=1.6, 
    #keysize=1.0, margin=c(25,20), cexRow=1.5, cexCol=1.6, 
    keysize=1.0, margin=c(25,20), cexRow=1.4, cexCol=2.0, 

    lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, .6, 0 ), lhei=c(.2, 3) , labRow = rowLabels_est)


    write.table(rowLabels_est, "rowLabels_est.txt", sep="\t")
    write.table(rowLabels_gene, "rowLabels_gene.txt", sep="\t")
    write.table(hm_gene_labels$rowInd, "hm_gene_labels_rowInd.txt", sep="\t")
    write.table(hm_est_labels$rowInd, "hm_est_labels_rowInd.txt", sep="\t")



   dev.off()


}


do_figure3 <- function() {
   datamatrix <<-read.table("annotated_data_matrix.txt", header=TRUE, sep="\t") 
   setwd("/dataset/bioinformatics_dev/active/qb_paper/qb_paper/figure_3") 
   jpeg(filename = "figure_3.jpg", width=800, height=800)
   plot(datamatrix[,"btau461_agbovine.csv"],datamatrix[,"umd3_agbovine.csv"],pch='.',cex=3.5,
            xlab="Bovine Genome (Btau 4.6.1)", ylab="Bovine Genome (UMD3)", cex.axis=1.5, cex.lab=1.5)
   text(8.18,12.35,"o",cex=4,col="red")
   title("EST Sequence Spectra (including only two genomes)", cex.main=2.0)

   dev.off()
}


#annotate_data()

#build_clustered_data()
#build_plot_data()
#do_figure2()
do_figure3()
#save.image(file="kmeans300_final.Rdata")
