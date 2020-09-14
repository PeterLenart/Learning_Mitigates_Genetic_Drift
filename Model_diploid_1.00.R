# working directory
setwd("c:/WORK/Research/Genetic drift and predation/")

# parameters
no_individuals <- 1000 # population size
init_no_alleles <- 100 # initial number of unique alleles
maxtime <- 100 # simulation time
no_remove <- 800 # number of individuals (prey) consumed each timestep
init_escape_prob <- 0.5 # initial probability of escaping predator attack
max_escape_prob <- 0.99 # maximum probability of escaping predator attack
rate_escape_prob <- 0.75 # proportional increase in escape probability (for learning; 0 for no learning)
no_simulations <- 10 # number of simulation runs

### with learning ###

allelesbeg <- array(0, c(2*no_individuals, no_simulations))
allelesend <- array(0, c(2*no_individuals, no_simulations))
no_diff_alleles <- array(0, c(no_simulations, maxtime))
totalremoved <- vector(mode="double", length=no_simulations)

for (s in 1:no_simulations) {
  
  escape <- array(init_escape_prob, c(no_individuals, maxtime))
  
  alleles <- array(0, c(no_individuals, 2, maxtime))
  alleles[,1,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  alleles[,2,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  allelesbeg[,s] <- c(alleles[,1,1], alleles[,2,1])
  
  no_diff_alleles[s,1] <- length(unique(c(alleles[,1,1], alleles[,2,1])))
  
  for (t in 2:maxtime) {
    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    
    actualprey <- 1:no_individuals
    repeat {
      attack <- sample(actualprey, 1) # random choice of prey
      if (runif(1) >= escape[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack, removed)
        actualprey <- setdiff(actualprey, attack)
      } else { # unsuccessful predator attack: increase in prey escape probability with learning
        escape[attack,t] <- min(max_escape_prob,(1+rate_escape_prob)*escape[attack,t-1])
      }
      if (successful_attacks >= no_remove) { break }
    }
    #print(length(removed))
    totalremoved <- totalremoved + length(removed)
    
    # determining indices of not removed individuals
    staying <- setdiff(1:no_individuals, removed)
    #print(length(staying))
    
    # staying individuals reproduce and produce no_remove offspring, to keep population size constant
    dalleles <- array(0, c(2, no_remove))
    for (i in 1:no_remove) {
      mates <- sample(staying, 2, replace=F)
      if (runif(1) < 0.5) {
        dalleles[1,i] <- alleles[mates[1],1,t-1]
      }
      else {
        dalleles[1,i] <- alleles[mates[2],1,t-1]
      }        
      if (runif(1) < 0.5) {
        dalleles[2,i] <- alleles[mates[1],2,t-1]
      }
      else {
        dalleles[2,i] <- alleles[mates[2],2,t-1]
      }        
    }
    alleles[,,t] <- alleles[,,t-1]
    alleles[removed,1,t] <- dalleles[1,]
    alleles[removed,2,t] <- dalleles[2,]
    
    escape[removed,t] <- rep(init_escape_prob, no_remove) # new individuals have no escape experience
    
    # calculate number of unique alleles
    no_diff_alleles[s,t] <- length(unique(c(alleles[,1,t], alleles[,2,t])))

  } # loop over t
  allelesend[,s] <- c(alleles[,1,maxtime], alleles[,2,maxtime])
  
} # loop over s

### without learning ###

allelesbeg2 <- array(0, c(2*no_individuals, no_simulations))
allelesend2 <- array(0, c(2*no_individuals, no_simulations))
no_diff_alleles2 <- array(0, c(no_simulations, maxtime))
totalremoved2 <- vector(mode="double", length=no_simulations)

for (s in 1:no_simulations) {
  
  escape <- array(init_escape_prob, c(no_individuals, maxtime))
  
  alleles <- array(0, c(no_individuals, 2, maxtime))
  alleles[,1,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  alleles[,2,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  allelesbeg2[,s] <- c(alleles[,1,1], alleles[,2,1])
  
  no_diff_alleles2[s,1] <- length(unique(c(alleles[,1,1], alleles[,2,1])))

  for (t in 2:maxtime) {
    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    
    actualprey <- 1:no_individuals
    repeat {
      attack <- sample(actualprey, 1) # random choice of prey
      if (runif(1) >= escape[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack, removed)
        actualprey <- setdiff(actualprey, attack)
      } 
      if (successful_attacks >= no_remove) { break }
    }
    #print(length(removed))
    totalremoved2 <- totalremoved2 + length(removed)
    
    # determining indices of not removed individuals
    staying <- setdiff(1:no_individuals, removed)
    #print(length(staying))
    
    # staying individuals reproduce and produce no_remove offspring, to keep population size constant
    dalleles <- array(0, c(2, no_remove))
    for (i in 1:no_remove) {
      mates <- sample(staying, 2, replace=F)
      if (runif(1) < 0.5) {
        dalleles[1,i] <- alleles[mates[1],1,t-1]
      }
      else {
        dalleles[1,i] <- alleles[mates[2],1,t-1]
      }        
      if (runif(1) < 0.5) {
        dalleles[2,i] <- alleles[mates[1],2,t-1]
      }
      else {
        dalleles[2,i] <- alleles[mates[2],2,t-1]
      }        
    }
    alleles[,,t] <- alleles[,,t-1]
    alleles[removed,1,t] <- dalleles[1,]
    alleles[removed,2,t] <- dalleles[2,]

    # calculate number of unique alleles
    no_diff_alleles2[s,t] <- length(unique(c(alleles[,1,t], alleles[,2,t])))

  } # loop over t
  allelesend2[,s] <- c(alleles[,1,maxtime], alleles[,2,maxtime])

} # loop over s

library("RColorBrewer")
col = brewer.pal(n = 3, name = "Paired")

ccex <- 1.5
ccex2 <- 1.5

par(mfrow=c(1,1))
xmax <- maxtime
ymax <- init_no_alleles
plot(1:maxtime, no_diff_alleles[1,], type = 'l', lwd = 3, col = col[1], xlab = 'Time', ylab = 'Number of different alleles', cex.axis = ccex, cex.lab = ccex, xlim = c(1, xmax), ylim = c(0, ymax))
lines(1:maxtime, no_diff_alleles2[1,], type = 'l', lwd = 3, col = col[2])
legend('bottomleft', c('with learning', 'without learning'), lty = 1, pch = -1, bty = 'o', lwd = 3, cex = ccex2, col = col)
for (s in 2:no_simulations) {
  lines(1:maxtime, no_diff_alleles[s,], type = 'l', lwd = 3, col = col[1])
  lines(1:maxtime, no_diff_alleles2[s,], type = 'l', lwd = 3, col = col[2])
}

# #png("hist.png", width=800, height=600)
# par(mfrow=c(2,2))
# ss <- 10
# hist(allelesbeg[,ss], main="With learning: t=0", xlab="Allele number")
# hist(allelesend[,ss], main="With learning: t=100", xlab="Allele number")
# hist(allelesbeg2[,ss], main="Without learning: t=0", xlab="Allele number")
# hist(allelesend2[,ss], main="Without learning: t=100", xlab="Allele number")
# par(mfrow=c(1,1))
# #dev.off()

print(totalremoved)
print(totalremoved2)
