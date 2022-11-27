# working directory
setwd("c:/WORK/Research/Genetic drift and predation/")

############
# parameters
############
no_individuals <- 1000 # population size
init_no_alleles <- 100 # initial number of unique alleles
maxtime <- 101 # simulation time
no_remove <- 300 # number of individuals (prey) consumed each timestep
init_escape_prob <- 0.5 # initial probability of escaping predator attack
max_escape_prob <- 0.99 # maximum probability of escaping predator attack
rate_escape_prob <- 0.75 # proportional increase in escape probability (for learning; 0 for no learning)
no_simulations <- 20 # number of simulation runs

### with learning ###

#allelesbeg <- array(0, c(2*no_individuals, no_simulations))
#allelesend <- array(0, c(2*no_individuals, no_simulations))
no_diff_alleles <- array(0, c(no_simulations, maxtime))

for (s in 1:no_simulations) {
  show(s)
  
  alleles <- array(0, c(no_individuals, 2, maxtime))
  alleles[,1,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  alleles[,2,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  #allelesbeg[,s] <- c(alleles[,1,1], alleles[,2,1])
  no_diff_alleles[s,1] <- length(unique(c(alleles[,1,1], alleles[,2,1])))
  escape_prob <- array(init_escape_prob, c(no_individuals, maxtime))
  
  for (t in 2:maxtime) {

    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    escape_prob_temp <- escape_prob[,t-1]
    
    actualprey <- 1:no_individuals
    repeat {
      attack <- sample(actualprey, 1) # random choice of prey
      if (runif(1) >= escape_prob[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack, removed)
        actualprey <- setdiff(actualprey, attack)
      } else { # unsuccessful predator attack: increase in prey escape probability with learning
        escape_prob_temp[attack] <- min(max_escape_prob,(1+rate_escape_prob)*escape_prob[attack,t-1])
        #escape_prob_temp[attack] <- escape_prob[attack,t-1] + rate_escape_prob*(max_escape_prob-escape_prob[attack,t-1]) 
      }
      if (successful_attacks >= no_remove) { break }
    }
    no_removed <- length(removed)
    no_staying <- length(actualprey)
    
    # no_staying individuals in actualprey reproduce and produce no_remove offspring, to keep population size constant
    alleles[,,t] <- alleles[,,t-1]
    dalleles <- array(0, c(2, no_remove))
    for (i in 1:no_removed) {
      mates <- sample(actualprey, 2, replace=F)
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
      jj <- removed[i]
      alleles[jj,1,t] <- dalleles[1,i] # inheritance
      alleles[jj,2,t] <- dalleles[2,i] # inheritance
      escape_prob_temp[jj] <- init_escape_prob
    }

    escape_prob[,t] <- escape_prob_temp 
    #escape_prob[removed,t] <- rep(init_escape_prob, no_removed) # new individuals have no escape experience

    # calculate number of unique alleles
    no_diff_alleles[s,t] <- length(unique(c(alleles[,1,t], alleles[,2,t])))

  } # loop over t
  #allelesend[,s] <- c(alleles[,1,maxtime], alleles[,2,maxtime])
  
} # loop over s

### without learning ###

#allelesbeg2 <- array(0, c(2*no_individuals, no_simulations))
#allelesend2 <- array(0, c(2*no_individuals, no_simulations))
no_diff_alleles2 <- array(0, c(no_simulations, maxtime))

for (s in 1:no_simulations) {
  show(s)
  
  alleles2 <- array(0, c(no_individuals, 2, maxtime))
  alleles2[,1,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  alleles2[,2,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  #allelesbeg2[,s] <- c(alleles[,1,1], alleles[,2,1])
  no_diff_alleles2[s,1] <- length(unique(c(alleles2[,1,1], alleles2[,2,1])))
  escape_prob2 <- array(init_escape_prob, c(no_individuals, maxtime))
  
  for (t in 2:maxtime) {
    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed2 <- NULL
    escape_prob_temp2 <- escape_prob2[,t-1]
    
    actualprey2 <- 1:no_individuals
    repeat {
      attack <- sample(actualprey2, 1) # random choice of prey
      if (runif(1) >= escape_prob2[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed2 <- c(attack, removed2)
        actualprey2 <- setdiff(actualprey2, attack)
      } 
      if (successful_attacks >= no_remove) { break }
    }
    no_removed <- length(removed2)
    no_staying <- length(actualprey2)
    
    # no_staying actualptey2 individuals reproduce and produce no_remove offspring, to keep population size constant
    alleles2[,,t] <- alleles2[,,t-1]
    dalleles <- array(0, c(2, no_remove))
    for (i in 1:no_removed) {
      mates <- sample(actualprey2, 2, replace=F)
      if (runif(1) < 0.5) {
        dalleles[1,i] <- alleles2[mates[1],1,t-1]
      }
      else {
        dalleles[1,i] <- alleles2[mates[2],1,t-1]
      }        
      if (runif(1) < 0.5) {
        dalleles[2,i] <- alleles2[mates[1],2,t-1]
      }
      else {
        dalleles[2,i] <- alleles2[mates[2],2,t-1]
      }        
      jj <- removed2[i]
      alleles2[jj,1,t] <- dalleles[1,i] # inheritance
      alleles2[jj,2,t] <- dalleles[2,i] # inheritance
      escape_prob_temp2[jj] <- init_escape_prob
    }

    # calculate number of unique alleles
    no_diff_alleles2[s,t] <- length(unique(c(alleles2[,1,t], alleles2[,2,t])))

  } # loop over t
  #allelesend2[,s] <- c(alleles[,1,maxtime], alleles[,2,maxtime])

} # loop over s

library("RColorBrewer")
col = brewer.pal(n = 3, name = "Paired")

ccex <- 1.5
ccex2 <- 1.5

suffix <- "_300rem_075esc_099max_default_20sim_noaging_accelerating";
sfile1 <- paste("diploid",suffix,".png",sep="")

png(file=sfile1,width=7.5,height=6,units="in",res=600)
par(mar=c(5,5,2,2))
xmax <- maxtime
ymax <- init_no_alleles
plot(1:maxtime, no_diff_alleles[1,], type = 'l', lwd = 3, col = col[1], xlab = 'Time', ylab = 'Number of different alleles', cex.axis = ccex, cex.lab = ccex, xlim = c(1, xmax), ylim = c(0, ymax))
lines(1:maxtime, no_diff_alleles2[1,], type = 'l', lwd = 3, col = col[2])
legend('bottomleft', c('with learning', 'without learning'), lty = 1, pch = -1, bty = 'o', lwd = 3, cex = ccex2, col = col)
for (s in 2:no_simulations) {
  lines(1:maxtime, no_diff_alleles[s,], type = 'l', lwd = 3, col = col[1])
  lines(1:maxtime, no_diff_alleles2[s,], type = 'l', lwd = 3, col = col[2])
}
text(60,3,'removed prey = 300', adj = c(0,0), cex = 1.5)

dev.off()

# #png("hist.png", width=800, height=600)
# par(mfrow=c(2,2))
# ss <- 10
# hist(allelesbeg[,ss], main="With learning: t=0", xlab="Allele number")
# hist(allelesend[,ss], main="With learning: t=100", xlab="Allele number")
# hist(allelesbeg2[,ss], main="Without learning: t=0", xlab="Allele number")
# hist(allelesend2[,ss], main="Without learning: t=100", xlab="Allele number")
# par(mfrow=c(1,1))
# #dev.off()

