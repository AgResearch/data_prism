###################################################################################
#  R scripts to load a dataset consisting of an alignment statistic (% identity) from
#  blast results of a set of 1K prob seqences against a number of genomes, and plot the results.
#  For a dataset containing over a million probes, this takes several hours to
#  run. In order to support incremental updates of plots etc,
#  the "time-consuming" results are assigned to global variables,
#  and thus saved as part of the workspace. Plots can then be adjusted by loading
#  the works-space back in and commenting out the time consuming global
#  variable assignments. 
###################################################################################
library(Heatplus)
library(RColorBrewer)
library("gplots")

###################################################################################
# edit these variables - the values given are just examples                       #
###################################################################################
data_folder="/dataset/bioinformatics_dev/active/qb_paper/qb_paper/figure_4"     
output_folder="/dataset/bioinformatics_dev/active/qb_paper/qb_paper/figure_4"   

heatmap_image_file="hm.jpg"
heatmap_image_file_subset="hm_subset.jpg"

mds_image_file="md1.jpg"
mds_image_file_subset="md_subset.jpg"

number_of_heatmap_row_labels=100
number_of_column_labels=30
num_clusters=10000
#num_clusters=200

selected_genomes<-c("EFGW73", "EFGW75", "EFGW7F", "EFGW7G", "EFGXW9", "EFGXX6", "EFG10DZ","NZ_CP010029.fasta")
###################################################################################
setwd(data_folder)

get_test_data<-function() {

    print("loading data")
    datamatrix1 <- read.table("probe_results.test1.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix2 <- read.table("probe_results.test2.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")

    print("merging data")
    datamatrix <<- cbind(datamatrix1, datamatrix2)
    colnames(datamatrix) <<- sub("probes_","",colnames(datamatrix))
    colnames(datamatrix) <<- sub(".out.gz","",colnames(datamatrix))
    print("done merging data - size is ")
    print(object.size(datamatrix))
}

get_data<-function() {
    # each of the results files loaded consists of an array of pct identity,
    # each probe is a row each column is a genome. Each partial file 
    # contains results from several genomes.

    print("loading data")
    datamatrix1 <- read.table("probe_results.00001.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix2 <- read.table("probe_results.00002.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix3 <- read.table("probe_results.00003.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix4 <- read.table("probe_results.00004.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix5 <- read.table("probe_results.00005.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix6 <- read.table("probe_results.00006.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix7 <- read.table("probe_results.00007.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix8 <- read.table("probe_results.00008.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix9 <- read.table("probe_results.00009.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix10 <- read.table("probe_results.00010.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix11 <- read.table("probe_results.00011.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix12 <- read.table("probe_results.00012.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")
    datamatrix13 <- read.table("probe_results.00013.txt_pctidentity.txt", header=TRUE, row.names=1, sep="\t")

    print("merging data")
    datamatrix <<- cbind(datamatrix1, datamatrix2, datamatrix3, datamatrix4, datamatrix5,datamatrix6, datamatrix7, datamatrix8, datamatrix9, datamatrix10, datamatrix11, datamatrix12, datamatrix13)
    colnames(datamatrix) <<- sub("probes_","",colnames(datamatrix))
    colnames(datamatrix) <<- sub(".out.gz","",colnames(datamatrix))
    print("done merging data - size is ")
    print(object.size(datamatrix))
}

do_clustering<-function() {
    # because of redundancy in the set of probes, the rows are clustered, yielding
    # a dataset with 1 row = 1 "probeset" (i.e. center of a cluster), 1 column = 1 genome
    print("clustering...")
    clustering <<- kmeans(datamatrix, num_clusters, iter.max=1000)
    print("done clustering")
}

do_hm_latest<-function(plotfile) {
    # present the results as a heatmap , rows = probesets,
    # columns = genomes
    print("doing heatmaps...")


    #3382    membrane protein
    # 3480    penicillin-insensitive murein endopeptidase
    cluster_gene<-read.table("cluster_genes.txt", header=TRUE, row.names=1, sep="\t")
    print(rownames(plotdata))
    print(cluster_gene)
    annotated_plotdata <- plotdata
    rownames(annotated_plotdata) <- as.integer(rownames(annotated_plotdata))
    print(rownames(annotated_plotdata))
    
    annotated_plotdata <- cbind(annotated_plotdata, as.character(cluster_gene[rownames(annotated_plotdata),"gene"]))
    colnames(annotated_plotdata) <- c(colnames(plotdata),"gene")

    #print(cluster_gene[1,"gene"])
    #print(plotdata["1",]) 

    #print(annotated_plotdata["1",]) 



    write.table(annotated_plotdata, "annotated_plotdata.txt", sep="\t")




    row_label_interval=max(1, floor(nrow(plotdata)/number_of_heatmap_row_labels))  # 1=label every location 2=label every 2nd location  etc
    col_label_interval=max(1, floor(ncol(plotdata)/number_of_column_labels))  # 1=label every location 2=label every 2nd location  etc

    #cm<-brewer.pal(9,"BuPu") # sequential
    cm<-brewer.pal(11,"Spectral")

    # set up a vector which will index the labels that are to be blanked out so that
    # only every nth row is labelled,
    # the rest empty strings, n=row_label_interval.
    #rowLabels <- rownames(plotdata)                   # these are the cluster numbers 
    rowLabels <- annotated_plotdata[,"gene"]           # these are gene descriptions 
    rowBlankSelector <- sequence(length(rowLabels))
    rowBlankSelector <- subset(rowBlankSelector, rowBlankSelector %% row_label_interval != 0)
                           # e.g. will get (2,3, 5,6, 8,9, ..)
                           # so we will only label rows 1,4,7,10,13 etc)

    # set up a vector which will index the labels that are to be blanked out so that
    # only every nth col is labelled,
    # the rest empty strings, n=col_label_interval.
    colLabels <- colnames(plotdata)
    colBlankSelector <- sequence(length(colLabels))
    colBlankSelector <- subset(colBlankSelector, colBlankSelector %% col_label_interval != 0)
                           # e.g. will get (2,3, 5,6, 8,9, ..)
                           # so we will only label rows 1,4,7,10,13 etc)

    jpeg(filename = "hm_internal.jpg", width=1300, height=1200) # with dendrograms

    # run the heatmap, just to obtain the clustering index - not the final plot
    hm_internal<-heatmap.2(as.matrix(plotdata),  scale = "none", dendrogram = "col",
         Colv = TRUE,
         trace = "none", breaks = c(90,95,97,99,99.125,99.25,99.375,99.5,99.625,99.75,99.875,100),
         col = cm , key=TRUE, density.info="none",
         keysize=1.0, margin=c(11,20), cexRow=1.5, cexCol=1.5,
         lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.7, 1.7, .6 ), lhei=c(.5, 3) , labRow = rowLabels)


    dev.off()

    # edit the re-ordered vector of row labels, obtained from the heatmap object, so that only
    # every nth label on the final plot has a non-empty string
    # this is for the internal distance matrix
    indexSelector <- hm_internal$rowInd[length(hm_internal$rowInd):1]
    indexSelector <- indexSelector[rowBlankSelector]
    rowLabels[indexSelector] = rep('',length(indexSelector))


    # edit the re-ordered vector of col labels, obtained from the heatmap object, so that only
    # every nth label on the final plot has a non-empty string
    # this is for the internal distance matrix
    indexSelector <- hm_internal$colInd[length(hm_internal$colInd):1]
    indexSelector <- indexSelector[colBlankSelector]
    colLabels[indexSelector] = rep('',length(indexSelector))


    # now do the final plot with all labels
    # https://www.rdocumentation.org/packages/gplots/versions/3.0.1/topics/heatmap.2
    jpeg(filename = plotfile, width=1400, height=1500) # with dendrograms
    #jpeg(filename = "figure_4_cluster_labels.jpg", width=1400, height=1500) # with dendrograms

    hm<<-heatmap.2(as.matrix(plotdata),  scale = "none", dendrogram = "col",
           Colv = TRUE,
           trace = "none", breaks = c(90,95,97,99,99.125,99.25,99.375,99.5,99.625,99.75,99.875,100),
           #col = cm , key=TRUE, density.info="none",
           col = cm , key=TRUE, density.info="none",
           #keysize=1.0, margin=c(40,60), cexRow=1.3, cexCol=1.3,
           #keysize=1.0, margin=c(20,20), cexRow=1.3, cexCol=1.3,
           #keysize=1.0, margin=c(20,20), cexRow=1.4, cexCol=2.0,
           #keysize=1.0, margin=c(20,40), cexRow=1.4, cexCol=2.0,
           #keysize=1.0, margin=c(20,40), cexRow=1.6, cexCol=2.0,
           keysize=1.0, margin=c(20,50), cexRow=1.7, cexCol=2.0,

           #lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, .8, 0 ), lhei=c(.5, 3) , labRow = rowLabels, labCol=colLabels)
           #lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, .8, 0 ), lhei=c(.2, 3) , labRow = rowLabels, labCol=colLabels)
           #lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.1, .8, 0 ), lhei=c(.2, 3) , labRow = rowLabels, labCol=colLabels)
           lmat=rbind(  c(4,3,0 ), c(2, 1, 0) ), lwid=c(.2, 1.2, 0 ), lhei=c(.2, 3) , labRow = rowLabels, labCol=colLabels)
    dev.off()

    # output the cluster info generated
    orderedmat=clustering$centers[hm$rowInd[length(hm$rowInd):1],]
    orderedmat=orderedmat[,hm$colInd[1:length(hm$colInd)]]  # NB, the column labels come out reversed !
    write.table(orderedmat,file="heatmap.dat",row.names=TRUE,sep="\t")
    write.table(clustering$cluster, file="heatmap_probsets.dat", row.names=TRUE, sep="\t")
    print("done heatmaps")
}

################### do mds of distances between genomes #############
do_mds<-function() {
    # distances between genomes are obtained as euclidean distances between the
    # columns of the clustered data. The resulting distance matrix is
    # embedded in a 2-D visualisation using multi-dimensional scaling
    
    print("doing mds plot...")

    distances = dist(as.matrix(t(clustering$centers)))
    fit <<- cmdscale(distances,eig=TRUE, k=2)

    jpeg(filename = mds_image_file, 800,800)
    #smoothScatter(fit$points, cex=0.8)
    #smoothScatter(fit$points, cex=0.8)
    text(fit$points, labels = rownames(fit$points), pos = 1, cex=0.8)
    dev.off()

    jpeg(filename = mds_image_file_subset, 800,800)
    #work out which ones to label
    comparison_points = as.vector(integer())
    comparison_points_patterns = selected_genomes
    #comparison_points_patterns = c("EFG10DQ", "EFG10DR")
    for (comparison_pattern in comparison_points_patterns){
             comparison_points = c(comparison_points, grep(comparison_pattern,rownames(fit$points), ignore.case = TRUE) )
    }
    point_labels=rep('', length(fit$points))
    point_labels[comparison_points] = rownames(fit$points)[comparison_points]

    smoothScatter(fit$points, cex=0.8)
    text(fit$points, labels = point_labels, pos = 1, cex=0.8)
    dev.off()
    print("done mds plot")
}

do_hm_subset<-function() {
    # this method used to make a heatmap from a subset of the data
    # (i.e. only for selected genomes)
    setwd(data_folder)
    print("loading workspace...")
    load("plots.RData")
    print("subsetting data...")
    plotdata<<-clustering$centers[,selected_genomes]
    do_hm(heatmap_image_file_subset)
}


create_cluster <- function(pos_x, pos_y, cluster_num, cluster_count) {
   # this method is used to manually adjust the "labelling clusters" , used to
   # label / annotate the mds plot - pass in the x, y coords of the new cluster and
   # the new cluster number and target membership. This method will move the
   # closest "cluster_count" points in to this cluster , adjust cluster memberships, and
   # re-calculate all the cluster centers. This is done to improve the labelling
   # of the mds plot
 
   
   new_center = c(pos_x, pos_y)
   dist_from_new <- fit$points  # just to get something with the right rownames. (We won't use the 2nd column)

   for (row_num in sequence(nrow(fit$points))) {
      v_data = as.numeric(fit$points[row_num,])
      d = (new_center - v_data) %*% (new_center - v_data)
      dist_from_new[row_num,1] = d
   }
   sorted_dist_from_new <- dist_from_new[order(dist_from_new[,1]),]
   points_to_move<-head(sorted_dist_from_new,cluster_count)

   # extend the sizes object
   sizes <<- c(sizes, 0)

   for(rowname in rownames(points_to_move)) {
      sizes[clusters[rowname]] <<- sizes[clusters[rowname]] - 1
      sizes[cluster_num] <<- sizes[cluster_num] + 1 
      clusters[rowname] <<- cluster_num      
   }
  
   rnames <- c(rownames(centers), cluster_num)
   centers <<- rbind(centers, c(pos_x, pos_y))
   rownames(centers) <<- rnames

   # recalculate centers
   for(icluster in sequence(length(sizes))) {
      centers[icluster, 1] <<- mean(fit$points[clusters==icluster,1])
      centers[icluster, 2] <<- mean(fit$points[clusters==icluster,2])
   }
 
   print(sizes)
   print(clusters)
   print(centers)
   print(points_to_move)
  
   #junk <- fit$points[order(fit$points[,1],fit$points[,2]),]
   #print(junk)
}


draw_distances_plot <- function() {
   # this method plots the mds view of the genomes distance
   # matrix, including adding non-overlapping labels of the data,
   # by grouping the points into "labelling clusters".
   jpeg(filename = mds_image_file, 800,800)
   #plot(fit$points, cex=0.7)
   smoothScatter(fit$points, cex=0.7)

   # in order to label, cluster the fit points and label the clusters
   #label_clustering <<- kmeans(fit$points, 17, iter.max=500)
   #label_clustering2 or 3 is best so far

   # re-assign these as we are going to edit the clustering later
   centers <<- label_clustering$centers
   sizes <<- label_clustering$size
   clusters <<- label_clustering$cluster

   # edits of clusters 
   #EFGXX2             458.649202  -72.600344
   #EFG10F4            460.131812  -62.477676
   #EFG10F5            541.825365 -107.405680
   #EFG10F3            551.082283  -96.994801
   create_cluster(503, -84, 16, 4) 

   #EFG10F1            171.201013 -239.211924
   #EFG10F2            173.437043 -240.404847
   #EFG10F6            173.768949 -240.268282
   #EFG10DV            186.055358 -261.110949
   #EFG10F7            186.635152 -261.145191
   create_cluster(178, -248, 17, 5) 


   #EFG10DX            -97.858328  -61.632361
   create_cluster(-98, -62, 18, 2) 


   #EFG10DW           -331.073329 -178.658835
   create_cluster(-331, -179, 19, 1) 

   #EFG10DS           -841.736298 -161.261213
   create_cluster(-842, -161, 20, 1) 

   #EFGW7C            -670.798716  216.858913
   #EFGW7B            -668.354855  215.990039
   create_cluster(-670, 216, 21, 2) 
   
   
   
   # label each non-trivial cluster with cluster count and 
   # draw a box around it. 
   xrad = (max(fit$points[,1]) - min(fit$points[,1]))/40
   yrad = (max(fit$points[,2]) - min(fit$points[,2]))/40

   if (min(sizes) == max(sizes)) {
      base = 1
   }
   else {
      base = max(sizes) - min(sizes)
   }
   
   for (i in sequence(length(sizes))) {
      if ( sizes[i] > 0 ) {
         #text(label_clustering$centers[ i,1 ] , 
         #           label_clustering$centers[ i,2 ] , 
         #              labels=label_clustering$size[i], cex = 1.0)
         xradius = ( xrad + 3*xrad * (sizes[i] - min(sizes))/ base ) / 2
         yradius = ( yrad + 3*yrad * (sizes[i] - min(sizes))/ base ) / 2
         rect(centers[ i,1 ] - xradius, centers[ i,2 ] - yradius,
              centers[ i,1 ] + xradius, centers[ i,2 ] + yradius) 
      }
   }

   #for (i in sequence(nrow(fit$points))) {
   #   if (label_clustering$size[label_clustering$cluster[i]] < 4) {
   #      text(fit$points[i,1], fit$points[i,2], rownames(fit$points)[i], cex=0.8)
   #   }
   #}


   # in case we want to label each cluster with the name of the genomes whose profile
   # is closest to the center of each cluster 
   #closest_dists = rep(NA,nrow(label_clustering$centers))
   #closest_rownums = rep(NA,nrow(label_clustering$centers))

   #for (center_num in sequence(nrow(label_clustering$centers))) {
   #   v_center = as.numeric(label_clustering$centers[center_num,])
   #   for (row_num in sequence(nrow(fit$points))) {
   #      v_data = as.numeric(fit$points[row_num,])
   #      d = (v_center - v_data) %*% (v_center - v_data)
   #      if(is.na(closest_dists[center_num])) {
   #         closest_dists[center_num] = d
   #         closest_rownums[center_num] = row_num
   #      }
   #      else if( d < closest_dists[center_num] ) {
   #         closest_dists[center_num] = d
   #         closest_rownums[center_num] = row_num
   #      }
   #   }
   #}

   # draw the labels to the clustered data
   # if we want to label with a representative genome
   #rownames=rownames(fit$points)[closest_rownums]
   #rownames=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O")

   
   for (i in sequence(length(sizes))) {
      yradius = ( yrad + 3*yrad * (sizes[i] - min(sizes))/ base ) / 2
      text(centers[ i,1 ] - 0.5 * xradius, centers[ i,2 ] + 0.8 * yradius, i, pos=3, cex=1.2)
   }


   # write out the clusters and points
   write.table(clusters,file="clusters.dat",row.names=TRUE,sep="\t")
   write.table(fit$points,file="points.dat",row.names=TRUE,sep="\t")

   dev.off()
}


#get_test_data()
#get_data()
#do_clustering()
#do_hm(hm.jpg)
print("loading workspace...")
load("plots2.RData")
number_of_heatmap_row_labels=80
#do_mds()
plotdata<-clustering$centers
do_hm_latest("figure_4.jpg")

#do_hm_subset()

#save.image(file="plots.RData")



