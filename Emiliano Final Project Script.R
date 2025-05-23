# R script for ECOL596G 
# Complex Systems & Networks - SPRING 2025
# written by Emiliano Calvo Alcaniz

## SET YOUR WORKING DIRECTORY:
setwd("C:/Desktop/Networks_final")
# Set this to wherever you put this folder

## ----------------
# Installing packages (*Make sure you already set your working directory*) ---------------------
install.packages("sna")
install.packages("intergraph")
install.packages("igraph")
install.packages("viridis")
install.packages("scales")
install.packages("bipartite")
install.packages("networkD3")
install.packages("circlize")
install.packages("snowboot")


## ----------------
# Loading in packages  ------------------
library(sna)
library(intergraph)
library(igraph)
library(viridis)
library(scales)
library(bipartite)  
library(networkD3) 
library(circlize)
library(snowboot)

## ----------------
# CAVIAR dataset IMPORTATION & calculations -----------------
# Load in the caviar dataset
caviar <- read.csv("CAVIAR_FULL.csv", header = TRUE)
# Turn it into a matrix
caviarm <- as.matrix(caviar[,1:110])

# Turn into an igraph form
caviar_igraph <- graph_from_adjacency_matrix(caviarm, mode = "directed", diag = TRUE ) 
# Simplifies the graph so that there aren't 200 edges between certain nodes
simplecaviar_igraph <- simplify(caviar_igraph)  

# Find the degree of each node in order to use when coloring the nodes in a plot
cdegree <- igraph::degree(simplecaviar_igraph, mode = "all") 
# Flatten and convert to numeric (cause it can't be used to plot color otherwise)
cdegree_numeric <- as.numeric(unlist(cdegree))  
# scaling the degree in order to color nodes
cd_max <- max(cdegree_numeric)
dc_scaled_calculations <- cdegree/cd_max
dc_scaled <- as.numeric(ceiling(dc_scaled_calculations * 100))
caviar_color <- rev(mako(100))
caviar_deg_colorscale <- caviar_color[dc_scaled]

# Creates the aesthetic parameter 
par(mfrow=c(1,2))
# because the main plot I'm using here is going to have two plots 
# next to each other because I'm comparing the degree to the 
# communication weights

plot(simplecaviar_igraph 
           # remove the labels because there's lots of nodes
     ,     vertex.label = NA                        
           # makes the color scale that i used before dependent on the degree
           # of the nodes
     ,     vertex.color = caviar_deg_colorscale
           # Shrink the arrow sizes because there's a lot of them
     ,     edge.arrow.size = 0.5
           )
# Add a title -- the "\n" knocks down the text by a line because it was too
# high in the plot
title(main = paste("\n                       ", "\n ",
  "\n Drug Smuggler Interactions colored by Degree"), cex.main = 1.5)  


# Create an object that contains the number of communications each individual
# received. ie. the edge weight
comm_rec <- as.numeric(rowSums(caviarm))
# Makes the colorscale dependent on the edge weight (communications received)
# rev() is used because it makes the dark colors heavier
comm_c_max <- max(comm_rec)
comm_c_scaled_calc <- comm_rec/comm_c_max
comm_scaled <- as.numeric(ceiling(comm_c_scaled_calc * 100))
commcaviar_colorscale <- caviar_color[comm_scaled]

communication_caviarplot <- plot(simplecaviar_igraph 
           # remove the labels because there's a lot of nodes
     ,     vertex.label = NA
           # makes the color scale dependent on the weight of communications
     ,     vertex.color = commcaviar_colorscale
           # shrinks the arrow sizes because there's a lot of them
     ,     edge.arrow.size = 0.5
)
# Add a title -- the "\n" knocks down the text by a line because it was too
# high in the plot
title(main = paste("\n                       ", "\n",
  "\n Drug Smuggler Interactions colored",  
                   "\n by Number of Communications Received"), cex.main = 1.5)  

## ----------------
#### CAVIAR NETWORK MEASURES  ------------
## YOU CAN JUST RUN THIS ENTIRE SECTION AND THE CONSOLE 
## WILL PRINT OUT ALL THE MEASURES (SO LONG AS YOU HAVE 
## ALREADY LOADED IN THE MATRICES)

# DEGREE
# Calculate degree (redundant, since it was calculated earlier for a plot)
caviar_deg <- igraph::degree(simplecaviar_igraph, mode = "in")
# Make it numeric so it can be used to calculate the mean
num_caviar_deg <- as.numeric(caviar_deg)
# Calculate the mean degree of the network
caviar_avg_deg <- mean(caviar_deg)
# Calculate it and print it out
cat("Average total degree of the directed network:", caviar_avg_deg, "\n")

## CLOSENESS

# Calculate the closeness
caviar_closeness <- igraph::closeness(simplecaviar_igraph)
# Calculate the mean of the closeness for the network
caviar_avg_closeness <- mean(caviar_closeness)
# Calculate it and print it out
cat("Average Closeness of this network:", caviar_avg_closeness, "\n")

## BETWEENNESS
caviar_betweenness <- igraph::betweenness(simplecaviar_igraph)
# Calculate the mean of the betweenness for the network
caviar_avg_betweenness <- mean(caviar_betweenness)
# Calculate it and print it out
cat("Average Betweenness of this network:", caviar_avg_betweenness, "\n")

## CLUSTERING COEFFICIENT

# Calculate the clustering coefficient-- in r, you use transitivity()
caviar_cc <- transitivity(simplecaviar_igraph, type = "global")
# Calculate it and print it out
cat("Clustering coefficient of this network:", caviar_cc, "\n")

# Create a random network [sample_gnm()] that has the same number of nodes 
# [vcount()] and edges [ecount()]
randomcaviar_cc1a <- sample_gnm(vcount(simplecaviar_igraph), ecount(caviar_igraph)  
                             , directed = TRUE, loops = FALSE)
# Calculate the Clustering coefficient of the random network
random_caviar_cc1 <- transitivity(randomcaviar_cc1a, type = "global")
# Calculate it and print it out
cat("Clustering coefficient of a RANDOM network:", random_caviar_cc1, "\n")



## MOTIFS

# Calculate the motifs of the dataset for groups of 3 individuals
caviar_mot <- motifs(simplecaviar_igraph, size = 3)
# Set the aesthetic parameters & layout
par(mar=c(0,1,0,1), oma = c(1,2,1,1), xpd=TRUE)# bottom, left, top, right
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17, 0)
              , nrow = 2, ncol = 18, byrow = T)
       , heights=c(4,1)) 
# set the colors for this graph, matching the color scale used earlier
color_distributions <- mako(16)
# create a barplot showing the number of motifs out of each possible motif
barplot(caviar_mot
        # color it using the mako color scale
        , col = color_distributions
        # idk what this does but it's necessary
        , names.arg = seq(1, 16)
)
# set the parameters to show the aesthetics going underneath
par(mar=c(0,0.5,0,0))
for (i in 0:15) {
  # This command gives the graph number i that is possible with
  # 3 nodes. 
  caviar_motif_graph <- graph_from_isomorphism_class(3, i)
  # Now plotting that:
  plot(caviar_motif_graph
       , edge.arrow.size = 0.5
       , edge.color = alpha("grey27", 0.5)
       , edge.width = 2
       , vertex.label.color = color_distributions
       , vertex.label.cex = 1
  )
}

## ASSORTATIVITY

# calculate the assortativity of the degree
caviar_assortativity <- assortativity(simplecaviar_igraph, values = caviar_deg, directed = TRUE)
# Calculate it and print it out
cat("Assortativity of this network:", caviar_assortativity, "\n")

# calculate the degree of the random graph from earlier
random_caviar_deg1 <- igraph::degree(randomcaviar_cc1a, mode = "in")
# Calculate the assortativity of the random graph
random_caviar_assortativity <- assortativity(randomcaviar_cc1a, values = random_caviar_deg1, directed = TRUE)
# Calculate it and print it out
cat("Assortativity of a RANDOM network:", random_caviar_assortativity, "\n")


## MODULARITY

# Calculate the modularity of this dataset
communities_caviar <- cluster_edge_betweenness(simplecaviar_igraph)
caviar_modularity <- modularity(communities_caviar)
# Calculate it and print it out
cat("Modularity of this network:", caviar_modularity, "\n")

## DIAMETER

# Calculate the diameter
caviar_diameter <- diameter(simplecaviar_igraph, directed = TRUE)
# Calculate it and print it out
cat("Diameter of this network:", caviar_diameter, "\n")



## ----------------
# EXERCISE-SPECIFIC CODE FOR THE CAVIAR DATASET - NOT FOR THE FINAL --------

# caviar_am <- as.matrix(as_adjacency_matrix(caviar_igraph))
# chordDiagram(caviar_am
#             , annotationTrack = "grid"      # This one erases the tick numbers
#             , grid.col = caviar_colorscale
#)

Ccloseness <- closeness(caviar_igraph)
Cbetweenness <- betweenness(caviar_igraph)
caviarplot_betweenness <- plot(Cbetweenness, col="blue4", bg="green4", 
                             lwd = 2, pch = 21, xlab = "Individual", ylab = "Betweenness")


# 4 random plots
color_closeness <- viridis(100)

randomcaviar1a <- sample_gnm(vcount(simplecaviar_igraph), ecount(simplecaviar_igraph)  
                         , directed = FALSE, loops = FALSE)
randomcaviar2a <- sample_gnm(vcount(simplecaviar_igraph), ecount(simplecaviar_igraph)
                         , directed = FALSE, loops = FALSE)
randomcaviar3a <- sample_gnm(vcount(simplecaviar_igraph), ecount(simplecaviar_igraph)
                         , directed = FALSE, loops = FALSE)
randomcaviar4a <- sample_gnm(vcount(simplecaviar_igraph), ecount(simplecaviar_igraph)
                          , directed = FALSE, loops = FALSE)

# Now we have three new networks, all resampled from the same distribution of
# random networks with the specified number of edges and vertices.
par(mfrow=c(2,2), mar = c(1,1,2,1))
mx_close <- max(closeness(caviar_igraph)
                , closeness(randomcaviar1a)
                , closeness(randomcaviar2a)
                , closeness(randomcaviar3a)
                , closeness(randomcaviar4a)
                , na.rm = T)

plot(randomcaviar1a
     , main = paste("V:", vcount(randomcaviar1a), " E:", ecount(randomcaviar1a))
     , vertex.color = color_closeness[100*closeness(randomcaviar1a)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(randomcaviar2a
     , main = paste("V:", vcount(randomcaviar2a), " E:", ecount(randomcaviar2a))
     , vertex.color = color_closeness[100*closeness(randomcaviar2a)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(randomcaviar3a
     , main = paste("V:", vcount(randomcaviar3a), " E:", ecount(randomcaviar3a))
     , vertex.color = color_closeness[100*closeness(randomcaviar3a)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)
plot(randomcaviar4a
     , main = paste("V:", vcount(randomcaviar4a), " E:", ecount(randomcaviar4a))
     , vertex.color = color_closeness[100*closeness(randomcaviar4a)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)



plot(caviar_igraph
     , main = paste("Original Caviar Network V:", vcount(caviar_igraph), " E:", ecount(caviar_igraph))
     , vertex.color = color_closeness[100*closeness(caviar_igraph)/mx_close]
     , vertex.label.cex = 0.7
     , vertex.label.color = "black"
     , vertex.label.family = "Arial"
     , edge.arrow.size = 0
     , edge.width = 2
)


color_closeness <- viridis(100)
plot(caviar_igraph
     , main = "Closeness"
     , vertex.color = color_closeness[100*closeness(caviar_igraph)/max(closeness(caviar_igraph))]
     , vertex.label.cex = 0.5
     , vertex.label.color = "black"
#     , vertex.label.family = "Arial"
     , rescale = FALSE
     , edge.width = 0.005
     , layout = caviar_coordinates
)

# heatmap(caviarm, col = rocket(50)) # HAS A BIG WHITE HORIZONTAL BAR BUT IT RUNS
# plot(as.network(caviarm)) 
# View(caviarm)
# caviarnetwork <- as.network(caviarm)
# plot(as.network(caviar), col = "mediumturquoise")


#### Weight adjustment
# caviar_list <- as.list(caviarm)
# omit_zeros_list <- caviar_list[caviar_list != 0]
# comm_to_weight <- function(x) {return(x/337)}
# E(caviar_igraph)$X <- lapply(E(caviar_igraph)
#                             , comm_to_weight)

# for infographic # par(mfrow=c(1,2)) # plot two graphs next to each other
# for infographic # par(mar = c(0,0,,0)) # bottom, left, top, right
# for infographic par(bg = '#fff9f2')
# for infographic:
# plot(caviarnetwork
#     , vertex.col= "#f9e1cf"
#     , edge.col = "#9c634f"
#     , main = "Network of Drug Smugglers"
#     , cex.main = 2)

# dev.off()



## ----------------
# LJ Dataset IMPORTATION & calculations  ------------------
ljedges <- read.csv("ljshrunkedges.csv", header = TRUE)
# Increases the 0s so r doesn't freak out when dividing by 0
ljm <- as.matrix(ljedges) + 1

# Making an igraph
ljg <- igraph::graph_from_edgelist(ljm, directed = TRUE)
simpleljg <- simplify(ljg, remove.multiple = TRUE, remove.loops = TRUE)
par(mfrow = c(1,1))
# Heatmap
lj_am <- as.matrix(as_adjacency_matrix(simpleljg))
nonzero <- rowSums(lj_am) != 0
lj_am_clean <- lj_am[nonzero, nonzero]
heatmap(lj_am_clean, rowv = NA, col = rocket(3))


plot(simpleljg
     , layout = layout_with_kk
     , rescale = TRUE
     , vertex.size = 4
     , vertex.label = NA
     , vertex.color = "lightsalmon1"
     , vertex.frame.color = "slateblue4"
     , edge.arrow.size = 0.25
)

title(main = paste("Network of Livejournal Friendships"), cex.main = 3)  


# Shrunken further to make the chord diagram
ljg_shrunk <- delete_vertices(ljg, v = c(150:1001))

lj_am_shrunk <- as.matrix(as_adjacency_matrix(
  ljg_shrunk
))

rownames(lj_am_shrunk) <- paste("I", c(1:vcount(ljg_shrunk)))

colnames(lj_am_shrunk) <- paste("I", c(1:vcount(ljg_shrunk)))

# circos.clear()     # Sometimes the diagram bugs out

lj_colorscale <- rocket(vcount(ljg_shrunk))

chordDiagram(lj_am_shrunk
             , annotationTrack = "grid"      # This one erases the tick numbers
             , grid.col = lj_colorscale
                 )


## ----------------
#### LJ  NETWORK MEASURES   -------------------------
## YOU CAN JUST RUN THIS ENTIRE SECTION AND THE CONSOLE 
## WILL PRINT OUT ALL THE MEASURES (SO LONG AS YOU HAVE 
## ALREADY LOADED IN THE MATRICES)

# DEGREE
lj_deg <- igraph::degree(simpleljg, mode = "in")
lj_avg_deg <- mean(lj_deg)
cat("Average total degree of the network:", lj_avg_deg, "\n")

# CLOSENESS
lj_closeness <- igraph::closeness(simpleljg)
lj_avg_closeness <- mean(lj_closeness)
cat("Average Closeness of this network:", lj_avg_closeness, "\n")

# BETWEENNESS
lj_betweenness <- igraph::betweenness(simpleljg, directed = TRUE)
lj_avg_betweenness <- mean(lj_betweenness)
cat("Average Betweenness of this network:", lj_avg_betweenness, "\n")

# CLUSTERING COEFFICIENT
lj_cc <- transitivity(simpleljg, type = "global")
cat("Clustering coefficient of this network:", lj_cc, "\n")

# Compared to a random network
randomlj_cc1a <- sample_gnm(vcount(simpleljg), ecount(simpleljg)  
                                , directed = TRUE, loops = FALSE)
random_lj_cc1 <- transitivity(randomlj_cc1a, type = "global")
cat("Clustering coefficient of a RANDOM network:", random_lj_cc1, "\n")

randomlj_cc2a <- sample_gnm(vcount(simpleljg), ecount(simpleljg)  
                            , directed = TRUE, loops = FALSE)


# MOTIFS
lj_mot <- motifs(simpleljg, size = 3)
par(mar=c(0,1,0,1), oma = c(1,2,1,1), xpd=TRUE)# bottom, left, top, right
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17, 0)
              , nrow = 2, ncol = 18, byrow = T)
       , heights=c(4,1)) 
ljcolor_distributions <- "salmon"
barplot(lj_mot
        , col = ljcolor_distributions
        , names.arg = seq(1, 16)
)
par(mar=c(0,0.5,0,0))
for (i in 0:15) {
  # This command gives the graph number i that is possible with
  # 3 nodes. 
  lj_motif_graph <- graph_from_isomorphism_class(3, i)
  # Now plotting that:
  plot(lj_motif_graph
       , edge.arrow.size = 0.5
       , edge.color = alpha("grey27", 0.5)
       , edge.width = 2
       , vertex.label.color = ljcolor_distributions
       , vertex.label.cex = 1
  )
}

par(mfrow = c(1,1))

# ASSORTATIVITY
lj_assortativity <- assortativity(simpleljg, values = lj_deg, directed = TRUE)
cat("Assortativity of this network:", lj_assortativity, "\n")

random_lj_deg1 <- igraph::degree(randomlj_cc1a, mode = "in")
random_lj_assortativity <- assortativity(randomlj_cc1a, values = random_lj_deg1, directed = TRUE)
cat("Assortativity of a RANDOM network:", random_lj_assortativity, "\n")


# MODULARITY
communities_lj <- cluster_edge_betweenness(simpleljg)
lj_modularity <- modularity(communities_lj)
cat("Modularity of this network:", lj_modularity, "\n")

# DIAMETER
lj_diameter <- diameter(simpleljg, directed = TRUE)
cat("Diameter of this network:", lj_diameter, "\n")





## ----------------
# Drosophila dataset IMPORTATION & calculations --------------

# Load in the csv file
dedges <- read.csv("drosophila_edgelist_shrunk.csv", header = TRUE)
# Turn it into a matrix
dedgem <- as.matrix(dedges)
# Turns the node numbers into characters so it doesn't bug out when i make it a graph
dedgem <- apply(dedgem, 2, as.character)
# turn the edgelist into a graph
dg <- graph_from_edgelist(dedgem)
# turn the graph into an adjacency matrix
dam <- as_adjacency_matrix(graph_from_edgelist(dedgem))
# I was just testing it here to make sure it would still work
dg <- graph_from_adjacency_matrix(dam, mode = "directed")

# View(as.matrix(dam))

# Aesthetics, this just makes the text italicized. I took it from a 
# stack overflow thing so I'm unclear on what I can remove without breaking it
italics <- lapply(rownames(data), function(x) bquote(italic(.(x))))
dr_color <- turbo(12)
dr_deg <- igraph::degree(dg, mode = "all")
dr_palette <- dr_color[dr_deg]

# PLOTTING 

dr_plot <- plot(intergraph::asNetwork(dg)
     , vertex.col = dr_palette
     , main = expression('Network of '~italic(Drosophila)~' Neurons')
     , cex.main = 2)

legend_values <- sort(unique(dr_deg))
legend_colors <- dr_color[legend_values]

legend("right",
       legend = legend_values,
       col = legend_colors,
       pch = 19,
       pt.cex = 2,
       title = "Degree",
       bty = "n")


# RANDOM PLOT
random_dr <- sample_gnm(vcount(dg), ecount(dg)  
                            , directed = TRUE, loops = FALSE)
dr_random_deg <- igraph::degree(random_dr, mode = "all")
dr_random_max <- max(dr_random_deg)
dr_random_color <- turbo(7)
dr_random_palette <- dr_random_color[dr_random_deg]

dr_random_plot <- plot(intergraph::asNetwork(random_dr)
                , vertex.col = dr_random_palette
                , main = expression('Random Network of '~italic(Drosophila)~' Neurons')
                , cex.main = 2)

legend_random_values <- sort(unique(dr_deg))
legend_random_colors <- dr_random_color[legend_random_values]

legend("right",
       legend = legend_random_values,
       col = legend_colors,
       pch = 19,
       pt.cex = 2,
       title = "Degree",
       bty = "n")



# get like a million warnings for dumb stuff so this clears it
warnings()

## ----------------
# PCA -- NOT FOR THE FINAL  -----

drosophilaloc <- cmdscale(dam)
drosophilax <- -drosophilaloc[, 1]
drosophilay <- -drosophilaloc[, 2] 

drosophila_am_result <- dam[, apply(dam, 2, var) != 0]

drosophila_pca_result <- prcomp(drosophila_am_result
                                , scale. = TRUE
                                )
drosophila_loadings <- drosophila_pca_result$rotation[, 1:2]

scaling_factor <- 1.5 * max(
  abs(c(drosophilax, drosophilay))) / max(abs(drosophila_loadings))

drosophila_loadings_scaled <- drosophila_loadings * scaling_factor

plot(drosophilax, drosophilay
     , main = expression('PCA of '~italic(Drosophila)~' Neurons')
     , cex.main = 2.5
     , xlab = "PC1", ylab = "PC2", axes = TRUE
     , asp = 1
     , pch = 19
     , col="skyblue3"
)

# arrows(0.05, 0, drosophila_loadings_scaled[,1]/3, drosophila_loadings_scaled[,2]/5,
#       length = 0.1, col = "chartreuse4", lwd = 2)



## ----------------
#### DROSOPHILA NETWORK MEASURES -----------------
## YOU CAN JUST RUN THIS ENTIRE SECTION AND THE CONSOLE 
## WILL PRINT OUT ALL THE MEASURES (SO LONG AS YOU HAVE 
## ALREADY LOADED IN THE MATRICES)

 # DEGREE
 dr_deg <- igraph::degree(dg, mode = "in")
 dr_avg_deg <- mean(dr_deg)
 cat("Average total degree of the network:", dr_avg_deg, "\n")
 
 # CLOSENESS
 dr_closeness <- igraph::closeness(dg)
 dr_avg_closeness <- mean(dr_closeness)
 cat("Average Closeness of this network:", dr_avg_closeness, "\n")
 
 # BETWEENNESS
 dr_betweenness <- igraph::betweenness(dg, directed = TRUE)
 dr_avg_betweenness <- mean(dr_betweenness)
 cat("Average Betweenness of this network:", dr_avg_betweenness, "\n")
 
 # CLUSTERING COEFFICIENT
 dr_cc <- transitivity(dg, type = "global")
 cat("Clustering coefficient of this network:", dr_cc, "\n")
 
 triad_census(dg) 
 
 
 # Compared to a random network
 randomdr_cc1a <- sample_gnm(vcount(dg), ecount(dg)  
                             , directed = TRUE, loops = FALSE)
 random_dr_cc1 <- transitivity(randomdr_cc1a, type = "global")
 cat("Clustering coefficient of a RANDOM network:", random_dr_cc1, "\n")

 
 # MOTIFS
 dr_mot <- motifs(dg, size = 3)
 par(mar=c(0,1,0,1), oma = c(1,2,1,1), xpd=TRUE)# bottom, left, top, right
 layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17, 0)
               , nrow = 2, ncol = 18, byrow = T)
        , heights=c(4,1)) 
 drcolor_distributions <- "skyblue3"
 barplot(dr_mot
         , col = drcolor_distributions
         , names.arg = seq(1, 16)
 )
 par(mar=c(0,0.5,0,0))
 for (i in 0:15) {
   # This command gives the graph number i that is possible with
   # 3 nodes. 
   dr_motif_graph <- graph_from_isomorphism_class(3, i)
   # Now plotting that:
   plot(dr_motif_graph
        , edge.arrow.size = 0.5
        #, edge.arrow.width = 1
        , edge.color = alpha("grey27", 0.5)
        , edge.width = 2
#        , vertex.label.family = "Arial"
        , vertex.label.color = drcolor_distributions
        , vertex.label.cex = 1
   )
 }
 
 par(mfrow = c(1, 1))
 
 # ASSORTATIVITY
 dr_assortativity <- assortativity(dg, values = dr_deg, directed = TRUE)
 cat("Assortativity of this network:", dr_assortativity, "\n")
 
 random_dr_deg1 <- igraph::degree(randomdr_cc1a, mode = "in")
 random_dr_assortativity <- assortativity(randomdr_cc1a, values = random_dr_deg1, directed = TRUE)
 cat("Assortativity of a RANDOM network:", random_dr_assortativity, "\n")
 
 
 # MODULARITY
 communities_dr <- cluster_edge_betweenness(dg)
 dr_modularity <- modularity(communities_dr)
 cat("Modularity of this network:", dr_modularity, "\n")
 
 # DIAMETER
dr_diameter <- diameter(dg, directed = TRUE)
 cat("Diameter of this network:", dr_diameter, "\n")
 


## ----------------
# FOCUS NETWORK -- WORKER INTERACTIONS dataset & calculations ------------------
## ----------------
# RAW DATA (JUST FOR ME - THE ENTIRE SECTION IS COMMENTS) -----
#Colony5CircleAggnT  <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Raw_Data/tracking_raw/Colony5_Circle_Aggn.csv")
#Colony5CircleAggnT <- setNames(Colony5CircleAggnT, c("ID", "Frames", "X", "Y", "Orientation",
#                                                     "SizeWidth.px", "SizeLeng.px", "Speed.Px.s", "Interpolated", "HeadX", "HeadY"))
#Colony5TubeAggnT  <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Raw_Data/tracking_raw/Colony5_Tube_Aggn.csv")
#Colony5TubeAggnT <- setNames(Colony5TubeAggnT, c("ID", "Frames", "X", "Y", "Orientation",
#                                                 "SizeWidth.px", "SizeLeng.px", "Speed.Px.s", "Interpolated", "HeadX", "HeadY"))
#
# BASELINE ASSAY
#Colony5CirclePreT <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Raw_Data/tracking_raw/Colony5_Circle_Pre.csv")
#Colony5CirclePreT <- setNames(Colony5CirclePreT, c("ID", "Frames", "X", "Y", "Orientation",
#                                                   "SizeWidth.px", "SizeLeng.px", "Speed.Px.s", "Interpolated", "HeadX", "HeadY"))
#Colony5TubePreT <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Raw_Data/tracking_raw/Colony5_Tube_Pre.csv")
#Colony5TubePreT <- setNames(Colony5TubePreT, c("ID", "Frames", "X", "Y", "Orientation",
#                                               "SizeWidth.px", "SizeLeng.px", "Speed.Px.s", "Interpolated", "HeadX", "HeadY"))

## ----------------
# DERIVED DATA -- dataset IMPORTATION ------------------

# INVADER Tube
Colony5TubeAggnNRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5TubeAggnNRMatrix.csv", row.names = 1, header = TRUE)
Colony5TubeAggnRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5TubeAggnRMatrix.csv", row.names = 1, header = TRUE)
# Matrix addition of the two matrices
Colony5TubeAggnMatrix <- as.matrix(Colony5TubeAggnNRMatrix + Colony5TubeAggnRMatrix)

colony5tubeaggn_igraph <- graph_from_adjacency_matrix(
  Colony5TubeAggnMatrix,
  mode = "directed",
  diag = TRUE,
)
colony5tubeaggnconnected_igraph <- delete_vertices(colony5tubeaggn_igraph, "id_359")

# INVADER Circle 
Colony5CircleAggnNRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5CircleAggnNRMatrix.csv", row.names = 1, header = TRUE) 
Colony5CircleAggnRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5CircleAggnRMatrix.csv", row.names = 1, header = TRUE) 
# Matrix addition of the two matrices
Colony5CircleAggnMatrix <- as.matrix(Colony5CircleAggnNRMatrix + Colony5CircleAggnRMatrix)

colony5circleaggn_igraph <- graph_from_adjacency_matrix(
  Colony5CircleAggnMatrix,
  mode = "directed",
  diag = TRUE,
)


# BASELINE Tube
Colony5TubePreNRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5TubePreNRMatrix.csv", row.names = 1, header = TRUE) 
Colony5TubePreRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5TubePreRMatrix.csv", row.names = 1, header = TRUE) 

# Matrix addition of the two matrices
Colony5TubePreMatrix <- as.matrix(Colony5TubePreNRMatrix + Colony5TubePreRMatrix)


colony5tubebaseline_igraph <- graph_from_adjacency_matrix(
  Colony5TubePreMatrix,
  mode = "directed",
  diag = TRUE,
)


# BASELINE Circle
Colony5CirclePreNRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5CirclePreNRMatrix.csv", row.names = 1, header = TRUE) 
Colony5CirclePreRMatrix <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/gchism/NestArchAggn/Derived_Data/Matrices/Colony5CirclePreRMatrix.csv", row.names = 1, header = TRUE) 
# Matrix addition of the two matrices
Colony5CirclePreMatrix <- as.matrix(Colony5CirclePreNRMatrix + Colony5CirclePreRMatrix)

colony5circlebaseline_igraph <- graph_from_adjacency_matrix(
  Colony5CirclePreMatrix,
  mode = "directed",
  diag = TRUE,
)


## ----------------
# Focal NETWORK MEASURES ----------------
## YOU CAN JUST RUN THIS ENTIRE SECTION AND THE CONSOLE 
## WILL PRINT OUT ALL THE MEASURES (SO LONG AS YOU HAVE 
## ALREADY LOADED IN THE MATRICES)

## CIRCLE BASELINE
fncirclebg <- simplify(colony5circlebaseline_igraph)

# DEGREE
fncb_deg <- igraph::degree(fncirclebg, mode = "in")
fncb_avg_deg <- mean(fncb_deg)
cat("Average total degree of the network:", fncb_avg_deg, "\n")

# CLOSENESS
fncb_closeness <- igraph::closeness(fncirclebg)
fncb_avg_closeness <- mean(fncb_closeness)
cat("Average Closeness of this network:", fncb_avg_closeness, "\n")

# BETWEENNESS
fncb_betweenness <- igraph::betweenness(fncirclebg, directed = TRUE)
fncb_avg_betweenness <- mean(fncb_betweenness)
cat("Average Betweenness of this network:", fncb_avg_betweenness, "\n")

# CLUSTERING COEFFICIENT
fncb_cc <- transitivity(fncirclebg, type = "global")
cat("Clustering coefficient of this network:", fncb_cc, "\n")

# ASSORTATIVITY
fncb_assortativity <- assortativity(fncirclebg, values = fncb_deg, directed = TRUE)
cat("Assortativity of this network:", fncb_assortativity, "\n")

# MODULARITY
communities_fncb <- cluster_edge_betweenness(fncirclebg)
fncb_modularity <- modularity(communities_fncb)
cat("Modularity of this network:", fncb_modularity, "\n")

# DIAMETER
fncb_diameter <- diameter(fncirclebg, directed = TRUE)
cat("Diameter of this network:", fncb_diameter, "\n")



## CIRCLE INVADER

fncircleag <- simplify(colony5circleaggn_igraph)

# DEGREE
fnca_deg <- igraph::degree(fncircleag, mode = "in")
fnca_avg_deg <- mean(fnca_deg)
cat("Average total degree of the network:", fnca_avg_deg, "\n")

# CLOSENESS
fnca_closeness <- igraph::closeness(fncircleag)
fnca_avg_closeness <- mean(fnca_closeness)
cat("Average Closeness of this network:", fnca_avg_closeness, "\n")

# BETWEENNESS
fnca_betweenness <- igraph::betweenness(fncircleag, directed = TRUE)
fnca_avg_betweenness <- mean(fnca_betweenness)
cat("Average Betweenness of this network:", fnca_avg_betweenness, "\n")

# CLUSTERING COEFFICIENT
fnca_cc <- transitivity(fncircleag, type = "global")
cat("Clustering coefficient of this network:", fnca_cc, "\n")


# ASSORTATIVITY
fnca_assortativity <- assortativity(fncircleag, values = fnca_deg, directed = TRUE)
cat("Assortativity of this network:", fnca_assortativity, "\n")

# MODULARITY
communities_fnca <- cluster_edge_betweenness(fncircleag)
fnca_modularity <- modularity(communities_fnca)
cat("Modularity of this network:", fnca_modularity, "\n")

# DIAMETER
fnca_diameter <- diameter(fncircleag, directed = TRUE)
cat("Diameter of this network:", fnca_diameter, "\n")



## TUBE BASELINE

fntubeb <- simplify(colony5tubebaseline_igraph)

# DEGREE
fntb_deg <- igraph::degree(fntubeb, mode = "in")
fntb_avg_deg <- mean(fntb_deg)
cat("Average total degree of the network:", fntb_avg_deg, "\n")

# CLOSENESS
fntb_closeness <- igraph::closeness(fntubeb)
fntb_avg_closeness <- mean(fntb_closeness)
cat("Average Closeness of this network:", fntb_avg_closeness, "\n")

# BETWEENNESS
fntb_betweenness <- igraph::betweenness(fntubeb, directed = TRUE)
fntb_avg_betweenness <- mean(fntb_betweenness)
cat("Average Betweenness of this network:", fntb_avg_betweenness, "\n")

# CLUSTERING COEFFICIENT
fntb_cc <- transitivity(fntubeb, type = "global")
cat("Clustering coefficient of this network:", fntb_cc, "\n")


# ASSORTATIVITY
fntb_assortativity <- assortativity(fntubeb, values = fntb_deg, directed = TRUE)
cat("Assortativity of this network:", fntb_assortativity, "\n")

# MODULARITY
communities_fntb <- cluster_edge_betweenness(fntubeb)
fntb_modularity <- modularity(communities_fntb)
cat("Modularity of this network:", fntb_modularity, "\n")

# DIAMETER
fntb_diameter <- diameter(fntubeb, directed = TRUE)
cat("Diameter of this network:", fntb_diameter, "\n")



## TUBE INVADER

fnta <- simplify(colony5tubeaggn_igraph)

# DEGREE
fnta_deg <- igraph::degree(fnta, mode = "in")
fnta_avg_deg <- mean(fnta_deg)
cat("Average total degree of the network:", fnta_avg_deg, "\n")

# CLOSENESS
fnta_closeness <- igraph::closeness(colony5tubeaggnconnected_igraph)
fnta_avg_closeness <- mean(fnta_closeness)
cat("Average Closeness of this network:", fnta_avg_closeness, "\n")

# BETWEENNESS
fnta_betweenness <- igraph::betweenness(fnta, directed = TRUE)
fnta_avg_betweenness <- mean(fnta_betweenness)
cat("Average Betweenness of this network:", fnta_avg_betweenness, "\n")

# CLUSTERING COEFFICIENT
fnta_cc <- transitivity(fnta, type = "global")
cat("Clustering coefficient of this network:", fnta_cc, "\n")

# ASSORTATIVITY
fnta_assortativity <- assortativity(fnta, values = fnta_deg, directed = TRUE)
cat("Assortativity of this network:", fnta_assortativity, "\n")

# MODULARITY
communities_fnta <- cluster_edge_betweenness(fnta)
fnta_modularity <- modularity(communities_fnta)
cat("Modularity of this network:", fnta_modularity, "\n")

# DIAMETER
fnta_diameter <- diameter(fnta, directed = TRUE)
cat("Diameter of this network:", fnta_diameter, "\n")




## ----------------
# Various t-tests ----------------

## BASELINES
#DEGREE
t.test(fncb_deg, fntb_deg)
#CLOSENESS
t.test(fncb_closeness,fntb_closeness)
#BETWEENNESS
t.test(fncb_betweenness, fntb_betweenness)

## INVADER ASSAYS
#DEGREE
t.test(fnca_deg, fnta_deg)
#CLOSENESS
t.test(fnca_closeness,fnta_closeness)
#BETWEENNESS
t.test(fnca_betweenness, fnta_betweenness)



## ----------------
# Robustness curve Circle ----------------
par(mfrow = c(1, 2))

fncirclebg_start <- simplify(colony5circlebaseline_igraph)
red_net_circle <- fncirclebg_start
eff_circle <- vector(mode = "numeric", length = vcount(fncirclebg_start))
shortest_dist_circle <- vector(mode = "numeric", length = vcount(fncirclebg_start))
comp_size_circle <- vector(mode = "numeric", length = vcount(fncirclebg_start))

for(i in 1:(vcount(red_net_circle)-1)) {
  eff_circle[i] <- global_efficiency(red_net_circle)
  shortest_dist_circle[i] <- mean_distance(red_net_circle)
  comp_size_circle[i] <- max(components(red_net_circle)$csize)
  vertex_deleted <- sample(V(red_net_circle), 1)
  if (vertex_deleted==13) (circle_node <- i)
  red_net_circle <- delete_vertices(red_net_circle, v = vertex_deleted)
}

## We've collected the data we want - note the above may take a while for a
## big network, and you may want to just calculate one measure to being with
## in that case.
## Now we plot the robustness curves:
cols <- cividis(3)
plot((shortest_dist_circle/max(shortest_dist_circle, na.rm=T))
     , main = "Robustness Curves of Baseline Circular Nest"
     , xlab = "Number of nodes removed +1"
     , ylab = "Measure"
     , ylim = c(0, 1)
     , pch = 19
     , col = cols[1]
)

points((eff_circle/max(eff_circle, na.rm = T))
       , pch = 19
       , col = cols[2]
)
points((comp_size_circle/max(comp_size_circle, na.rm = T))
       , pch = 19
       , col = cols[3]
)
legend("bottomleft"
       , legend = c("Shortest distance", "Efficiency", 
                     "Largest component")
       , col = cols, pch = 19
)

# Robustness curve Tube ----------------

fntubebg_start <- simplify(colony5tubebaseline_igraph)
red_net_tube <- fntubebg_start
eff_tube <- vector(mode = "numeric", length = vcount(fntubebg_start))
shortest_dist_tube <- vector(mode = "numeric", length = vcount(fntubebg_start))
comp_size_tube <- vector(mode = "numeric", length = vcount(fntubebg_start))

for(i in 1:(vcount(red_net_tube)-1)) {
  eff_tube[i] <- global_efficiency(red_net_tube)
  shortest_dist_tube[i] <- mean_distance(red_net_tube)
  comp_size_tube[i] <- max(components(red_net_tube)$csize)
  vertex_deleted <- sample(V(red_net_tube), 1)
  if (vertex_deleted==13) (tube_node <- i)
  red_net_tube <- delete_vertices(red_net_tube, v = vertex_deleted)
}

## We've collected the data we want - note the above may take a while for a
## big network, and you may want to just calculate one measure to being with
## in that case.
## Now we plot the robustness curves:
cols <- cividis(3)
tube_robustness_plot <- plot((shortest_dist_tube/max(shortest_dist_tube, na.rm=T))
     , main = "Robustness Curves of Baseline Tube Nest"
     , xlab = "Number of nodes removed +1"
     , ylab = "Measure"
     , ylim = c(0, 1)
     , pch = 19
     , col = cols[1]
)

points((eff_tube/max(eff_tube, na.rm = T))
       , pch = 19
       , col = cols[2]
)
points((comp_size_tube/max(comp_size_tube, na.rm = T))
       , pch = 19
       , col = cols[3]
)
legend("bottomleft"
       , legend = c("Shortest distance", "Efficiency", 
                    "Largest component")
       , col = cols, pch = 19
)



## ----------------
# Focal Network Plotting, etc ---------------
par(mfrow = c(1,1))

# Aesthetics
fn_color <- cividis(100)
par(mfrow = c(2, 2))

## SETUP
# CIRCLE BASELINE
fncb_closeness <- igraph::closeness(fncirclebg)
# Flatten and convert to numeric (cause it can't be used to plot color otherwise)
cbc_numeric <- as.numeric(unlist(fncb_closeness))  
# scaling the degree in order to color nodes
cbc_max <- max(cbc_numeric)
cbc_scaled_calculations <- cbc_numeric/cbc_max
cbc_scaled <- as.numeric(ceiling(cbc_scaled_calculations * 100))
cbc_colorscale <- fn_color[cbc_scaled]

# CIRCLE INVADER
fnca_closeness <- igraph::closeness(colony5circleaggn_igraph)
# Flatten and convert to numeric (cause it can't be used to plot color otherwise)
cac_numeric <- as.numeric(unlist(fnca_closeness))  
# scaling the degree in order to color nodes
cac_max <- max(cac_numeric)
cac_scaled_calculations <- cac_numeric/cac_max
cac_scaled <- as.numeric(ceiling(cac_scaled_calculations * 100))
cac_colorscale <- fn_color[cac_scaled]

# TUBE BASELINE
fntb_closeness <- igraph::closeness(colony5tubebaseline_igraph)
# Flatten and convert to numeric (cause it can't be used to plot color otherwise)
tbc_numeric <- as.numeric(unlist(fntb_closeness))  
# scaling the degree in order to color nodes
tbc_max <- max(tbc_numeric)
tbc_scaled_calculations <- tbc_numeric/tbc_max
tbc_scaled <- as.numeric(ceiling(tbc_scaled_calculations * 100))
tbc_colorscale <- fn_color[tbc_scaled]

# TUBE INVADER
fnta_closeness <- igraph::closeness(colony5tubeaggnconnected_igraph)
# Flatten and convert to numeric (cause it can't be used to plot color otherwise)
tac_numeric <- as.numeric(unlist(fnta_closeness))  
# scaling the degree in order to color nodes
tac_max <- max(tac_numeric)
tac_scaled_calculations <- tac_numeric/tac_max
tac_scaled <- as.numeric(ceiling(tac_scaled_calculations * 100))
tac_colorscale <- fn_color[tac_scaled]


# Plotting
# Circle Baseline Plot
par(mar=c(.5, 0.2, 2, 0.2), mfrow=c(2,2))

plot(simplify(colony5circlebaseline_igraph)
     , main = "(a) Circle Baseline Assay"
     , rescale = TRUE
     , vertex.label = NA
     , vertex.color = cbc_colorscale
     , edge.arrow.size = 0.5)

# Circle Invader Plot
plot(simplify(colony5circleaggn_igraph)
     , main = "(b) Circle Invader Assay"
     , rescale = TRUE
     , vertex.label = NA
     , vertex.color = cac_colorscale
     , edge.arrow.size = 0.5)

# Tube Baseline Plot
plot(simplify(colony5tubebaseline_igraph)
     , main = "(c) Tube Baseline Assay"
     , rescale = TRUE
     , vertex.label = NA
     , vertex.color = tbc_colorscale
     , edge.arrow.size = 0.5)


# Tube Invader Plot
plot(simplify(colony5tubeaggnconnected_igraph)
     , main = "(d) Tube Invader Assay"
     , rescale = TRUE
     , vertex.label = NA
     , vertex.color = tac_colorscale
     , edge.arrow.size = 0.5)

mtext("Interaction Networks \n Colored by Closeness"
      , cex = 2
      , side = 3
      , line = - 29
      , outer = TRUE)


# END ---------

