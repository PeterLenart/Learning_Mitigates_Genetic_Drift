# working directory
setwd("c:/Users/lberec/source/repos/LearningGeneticDrift_IBM/bin/Debug/netcoreapp3.1/")

suffix <- "_immediate_10_noaging_accelerating.txt"
file1 <- paste("heatmap",suffix,sep="")
file5 <- paste("learn_dif",suffix,sep="")
file6 <- paste("GT_dif",suffix,sep="")
file7 <- paste("GTA_dif",suffix,sep="")
file8 <- paste("AD_dif",suffix,sep="")

grid_x = 21 
grid_y = 21
start_rate = 0
end_rate = 0.9
start_remove = 0
end_remove = 900

heat_data <- read.table(file8, header = FALSE, sep = " ", dec = ".")
heat_data <- as.matrix(heat_data[,1:grid_y])
heat_data <- t(heat_data)
#heat_data[1,] <- 0                       


library(lattice)

fun.1 = function(x){x*((end_rate - start_rate)/(grid_y-1))} # collumn names
fun.2 = function(x){x*((end_remove - start_remove)/(grid_x-1))} # row names

collumns <- data.frame(x = c(start_rate:(grid_y-1)), y = NA)
collumns$y <- fun.1(collumns$x)
rows <- data.frame(x = c(start_remove:(grid_x-1)), y = NA)
rows$y <- fun.2(rows$x)

colnames(heat_data) <- collumns$y
rownames(heat_data) <- rows$y 

levelplot(heat_data, 
          
          #main="Difference in number of alleles after 50 generations (default)", 
          main="Difference in generation time (default)", 
          xlab=list(label="Number of prey removed per time step", cex=1.5),
          ylab=list(label="Rate of change of escape probability", cex=1.5),
          col.regions = terrain.colors(35), scales=list(x=list(cex=1.3, rot=90), y=list(cex=1.3)),
          colorkey=list(labels=list(cex=1, font=2)))

# run before new simulation

#rm(list = ls())
