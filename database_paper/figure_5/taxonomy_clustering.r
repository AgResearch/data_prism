# 
setwd("/dataset/hiseq/scratch/postprocessing/")

get_command_args <- function() {
   args=(commandArgs(TRUE))
   if(length(args)!=1 ){
      #quit with error message if wrong number of args supplied
      print('Usage example : Rscript --vanilla  taxonomy_clustering.r run_name=160623_D00390_0257_AC9B0MANXX')
      print('args received were : ')
      for (e in args) {
         print(e)
      }
      q()
   }else{
      print("Using...")
      # seperate and parse command-line args
      for (e in args) {
         print(e)
         ta <- strsplit(e,"=",fixed=TRUE)
         switch(ta[[1]][1],
            "run_name" = run_name <- ta[[1]][2]
         )
      }
   }
   return(run_name)
}

get_clusters <- function(datamatrix) {
   num_clust=20   # all 
   num_clust=20   # euk only
   clustering <- kmeans(datamatrix, num_clust, iter.max=500)
   tcenters=t(clustering$centers)
   distances <- dist(tcenters)
   sample_species_table<-read.table("sample_species.txt", header=TRUE, row.names=1, sep="\t")
   species_config_table<-read.table("species_config.txt", header=TRUE, row.names=1, sep="\t")
   sample_names <- strsplit(rownames(t(datamatrix)), split="_")
   get_name <- function(split_result) unlist(split_result)[5]
   sample_names <- sapply(sample_names, get_name) 
  
   # look up the species in the sample-species table , using sample name as key
   get_species <- function(sample_name) as.vector(sample_species_table[sample_name, "species"])[1]
   species_names <- sapply(sample_names, get_species)
   species_names <- sapply(species_names, tolower)

   # for each species in the species_config table, use the regexp 
   # in that table to locate the matching samples, and hence set the 
   # appropriate plot colour 
   point_colours <- rep("black", length(species_names)) 
   point_colours[1:(length(point_colours)-8)]<-NA

   point_symbols <- rep('?', length(species_names))
   point_count <- 1
   for(rowname in rownames(species_config_table)) {
      regexp <-  toString(species_config_table[rowname, "regexp"])
      colour <- paste("#", toString(species_config_table[rowname, "colour"]),sep="")
      symbol <- toString(species_config_table[rowname, "symbol"]) 
      point_colours[ grep(regexp,species_names, ignore.case = TRUE) ] <- colour
      point_symbols[ grep(regexp,species_names, ignore.case = TRUE) ] <- symbol
   }

   distances = dist(as.matrix(t(datamatrix)))

   fit <- cmdscale(distances,eig=TRUE, k=2)

   results=list()
   results$fit = fit
   results$point_symbols = point_symbols
   results$point_colours = point_colours
   results$sample_names = sample_names

   return(results)
}


stroverlap <- function(x1,y1,s1, x2,y2,s2) {
   # ref : https://stackoverflow.com/questions/6234335/finding-the-bounding-box-of-plotted-text
   # example : stroverlap(.5,.5,"word", .6,.5, "word")
   sh1 <- strheight(s1)
   sw1 <- strwidth(s1)
   sh2 <- strheight(s2)
   sw2 <- strwidth(s2)

   overlap <- FALSE
   if (x1<x2) 
     overlap <- x1 + sw1 > x2
   else
     overlap <- x2 + sw2 > x1

   if (y1<y2) 
     overlap <- overlap && (y1 +sh1>y2)
   else
     overlap <- overlap && (y2+sh2>y1)

   return(overlap)
}

plot_data<-function(filename, plot_title, save_prefix, cex_label, exclude_missing) {
   datamatrix<-read.table(filename, header=TRUE, row.names=1, sep="\t")
   clusters=get_clusters(datamatrix)

   if(exclude_missing) {
      fit_points = subset(clusters$fit$points, ! is.na(clusters$point_colours))
      point_symbols = subset(clusters$point_symbols, ! is.na(clusters$point_colours))
      sample_names = subset(clusters$sample_names, ! is.na(clusters$point_colours))
      point_colours = subset(clusters$point_colours, ! is.na(clusters$point_colours))
   }
   else {
      fit_points = clusters$fit$points 
      point_symbols = clusters$point_symbols
      sample_names = clusters$sample_names
      point_colours = clusters$point_colours
   }

   #print("symbols and colours")
   #print(point_symbols)
   #print(point_colours)

   #plot.default(fit_points, col=point_colours, cex=1.5, pch=point_symbols)
   #plot.default(fit_points, col=point_colours, cex=1.5, pch=point_symbols, xlab="", ylab="", cex.axis=1.2, cex.lab=1.2)
   #plot.default(fit_points, col=point_colours, cex=1.5, pch=point_symbols, xlab="", ylab="", cex.axis=1.2, cex.lab=1.2)
   plot.default(fit_points, col=point_colours, cex=1.5, pch=point_symbols, xlab="", ylab="", cex.axis=1.5, cex.lab=1.5)

   sample_names[1:(length(sample_names)-8)]<-NA
   title(plot_title, cex.main=1.5)

   # this next block of code is to do with avoiding over-plotting the labels
   #text(clusters$fit$points, labels = clusters$sample_names, pos = 1, cex=1.5)
   #print(clusters$fit$points)
   #print(length( clusters$sample_names ))

   # we will look for label overlaps , and form the overlapping labels into groups, 
   # and label each group with just one of the labels. Then we will also
   # put up a key in the top left, with the definition of each group 
   # (There will be some placements of labels which will still result in over-plotting,
   # - can adjust the sensitivity/specificity of the overlap detection to handle this
   # (and also tolerate some over-plotting) )
   # this section defines the groups. . . . 
   plot_group_labels = vector("list", 0) 
   plot_group_pos = vector("list",0)
   for(i in 1:length( sample_names )) {
      if (! is.na(sample_names[i])) {
         if( length( plot_group_labels ) == 0 ) {
            plot_group_labels = c(plot_group_labels, sample_names[i])
            plot_group_pos = c(plot_group_pos,"") # need to be careful extending lists, to contain something that is itself a list(else you will flatten the new member)
            plot_group_pos[[length(plot_group_pos)]] = c( fit_points[i,1], fit_points[i,2] )
         }
         else {
            assigned_to_group = FALSE
            for(j in 1:length( plot_group_labels )) {
               # if this label overlaps a plot group , append the name to the group, and we assigned this label to a group 
               if(stroverlap( fit_points[i,1], fit_points[i,2] , sample_names[i],
                  plot_group_pos[[j]][1], plot_group_pos[[j]][2], plot_group_labels[[j]] )) {
                  plot_group_labels[[j]] = c( plot_group_labels[[j]] , sample_names[i]) 
                  assigned_to_group = TRUE
                  next   
               }
            } #for each plot group 
            if( ! assigned_to_group ) {
               # label did not overlap any groups so make a singleton for it
               plot_group_labels = c(plot_group_labels, sample_names[i])
               plot_group_pos = c(plot_group_pos,"") # need to be careful extending lists, to contain something that is itself a list(else you will flatten the new member)
               plot_group_pos[[length(plot_group_pos)]] = c( fit_points[i,1], fit_points[i,2] )
            } # not found 
         } # plot group has been initialiased
      } # if one of the labels we are to plot 
   } # for each point
   print(plot_group_labels)
   print(plot_group_pos)
                       

   #for(i in 1:length( clusters$sample_names )) {
   #   #print(clusters$fit$points[i])
   #   #print(clusters$sample_names[i])
   #   print(paste(clusters$fit$points[i,1], clusters$fit$points[i,2], clusters$sample_names[i]))
   #   text(clusters$fit$points[i,1], y=clusters$fit$points[i,2], labels = clusters$sample_names[i], pos = 4, cex=cex_label)
   #}

   # label each "label-group" , with just one of the labels
   for(i in 1:length( plot_group_labels )) {
      if (length( plot_group_labels[[i]] ) ==1 ) {
         text(plot_group_pos[[i]][1], y=plot_group_pos[[i]][2], labels = plot_group_labels[[i]][1], adj=c(0,.5), cex=cex_label)
      }
      else {
         text(plot_group_pos[[i]][1], y=plot_group_pos[[i]][2], labels = paste(plot_group_labels[[i]][1], "(etc.)", sep=" "), adj=c(0,.5), cex=cex_label)
      }
   }


   # emit a key, defining the label groups (if there are any)
   top_left=c( 1.02 * min( fit_points[,1] ) , .98 * max( fit_points[,2] ) )
   key_row_count = 0
   for(i in 1:length( plot_group_labels )) {
      if (length( plot_group_labels[[i]] ) > 1 ) {
         key_string = paste(plot_group_labels[[i]][1], "(etc.) also includes nearby samples:", sep=" ")
         key_string = paste(key_string,  paste( plot_group_labels[[i]][2:length( plot_group_labels[[i]])], collapse=","))
         text(top_left[1] , y=top_left[2] - key_row_count * strheight(key_string) * 1.7, labels=key_string, cex=cex_label, pos=4)
         key_row_count = key_row_count + 1
      }
   }



   write.table(fit_points,file=paste(save_prefix,run_name,".txt",sep=""),row.names=TRUE,sep="\t")
}


run_name<<-get_command_args()


#jpeg(filename = paste("taxonomy_clustering_",run_name,".jpg",sep=""), 800, 800) # this setting used for doing just one of the three plots below 


# these settings no longer used  - concatenate images outside this script not 
#jpeg(filename = paste("taxonomy_clustering_",run_name,".jpg",sep=""), 800, 2400)
#par(mfrow=c(3, 1))
#cex_label=1.5   


cex_label=1.5   # this setting used for doing just one of the plots
jpeg(filename = paste("euk_taxonomy_clustering_",run_name,".jpg",sep=""), 800, 800) # this setting used for doing just one of the three plots below 
plot_data("eukaryota_information.txt", "Clustering of Blast Eukaryota-Hit Profiles", "Clustering-of-Blast-Eukaryota-Hit-Profiles-", cex_label, TRUE)
dev.off()

#jpeg(filename = paste("all_taxonomy_clustering_",run_name,".jpg",sep=""), 800, 800) # this setting used for doing just one of the three plots below 
#plot_data("all_information.txt", "Clustering of Blast All-Hit Profiles", "Clustering-of-Blast-All-Hit-Profiles-", cex_label, TRUE)
#dev.off()

#jpeg(filename = paste("xno_taxonomy_clustering_",run_name,".jpg",sep=""), 800, 800) # this setting used for doing just one of the three plots below 
#plot_data("all_information_xnohit.txt", "Clustering of Blast All-Hit Profiles (Excluding 'no hit')", "Clustering-of-Blast-All-Hit-Profiles-Excluding-no-hit-", cex_label, TRUE)
#dev.off()



