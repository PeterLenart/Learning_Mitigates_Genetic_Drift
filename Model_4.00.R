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
no_simulations <- 10 # number of simulation runs
haploid_asexual <- T # haploid asexual (T) or haploid sexual (F) prey
aging <- F # finite lifespan set on (T) or off (F)
lifespan <- 10 # maximum age if aging = T
random_aging <- F # finite lifespan (F) or random mortality every year (T)
immediate_update <- F # original update once per time step (F) or after every unsuccessful attack
# age distribution
age_time_sample <- 50
age_limit <- 16
############

### with learning ###

#allelesbeg <- array(0, c(no_individuals, no_simulations))
#allelesend <- array(0, c(no_individuals, no_simulations))
no_diff_alleles <- array(0, c(no_simulations, maxtime))
#totalremoved <- rep(0, no_simulations)
#totalattacks <- rep(0, no_simulations)
age <- rep(0, no_individuals)
#totalremovedx <- rep(0, no_simulations)
age_at_death <- NULL
age_dist <- array(0, c(age_limit, no_simulations))

for (s in 1:no_simulations) {
  
  age <- sample(0:(lifespan-1), no_individuals, replace=T)
  alleles <- array(0, c(no_individuals, maxtime))
  alleles[,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  #allelesbeg[,s] <- alleles[,1]
  no_diff_alleles[s,1] <- length(unique(alleles[,1]))
  escape_prob <- array(init_escape_prob, c(no_individuals, maxtime))
  
  for (t in 2:maxtime) {
    # saving age distribution
    if (t == age_time_sample) {
      for (i in 1:(age_limit-1)) {
        age_dist[i,s] <- sum(age == i-1)
      }
      age_dist[age_limit,s] <- sum(age >= (age_limit-1))
      show(sum(age_dist))
    }
    
    # aging
    age <- age + 1
    
    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    escape_prob_temp <- escape_prob[,t-1]
    
    actualprey <- 1:no_individuals
    repeat {
      attack <- sample(actualprey, 1) # random choice of prey
      #totalattacks[s] <- totalattacks[s] + 1
      if (runif(1) >= escape_prob[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack, removed)
        actualprey <- setdiff(actualprey, attack)
        age_at_death <- c(age_at_death, age[attack])
      } else { # unsuccessful predator attack: increase in prey escape probability with learning
        if (immediate_update == F) { # between timestep update
          escape_prob_temp[attack] <- min(max_escape_prob,(1+rate_escape_prob)*escape_prob[attack,t-1])
          #escape_prob_temp[attack] <- escape_prob[attack,t-1] + rate_escape_prob*(max_escape_prob-escape_prob[attack,t-1]) 
        }
        else { # within timestep update
          escape_prob[attack,t-1] <- min(max_escape_prob,(1+rate_escape_prob)*escape_prob[attack,t-1]) 
          #escape_prob[attack,t-1] <- escape_prob[attack,t-1] + rate_escape_prob*(max_escape_prob-escape_prob[attack,t-1]) 
        }
      }
      if (successful_attacks >= no_remove) { break }
    }
    #show(successful_attacks)
    #show(length(removed))
    #totalremoved[s] <- totalremoved[s] + length(removed)
    
    # remove surviving old individuals
    if (aging) {
      apC <- length(actualprey)
      if(!random_aging) {
        for (i in apC:-1:1) {
          ii <- actualprey[i]
          if (age[ii] > lifespan) {
            age_at_death <- c(age_at_death, age[ii])
            removed <- c(ii, removed)
            actualprey <- setdiff(actualprey, ii)
          }
        }
      }
      else { # random_aging == T
        for (i in apC:-1:1) {
          ii <- actualprey[i]
          if (runif(1) < 1/lifespan) {
            age_at_death <- c(age_at_death, age[ii])
            removed <- c(ii, removed)
            actualprey <- setdiff(actualprey, ii)
          }
        }
      }
      #totalremovedx[s] <- totalremovedx[s] + length(removed2)
      #staying2 <- setdiff(staying, removed2) # finally staying
      #removed <- c(removed, removed2) # total removed
    }
    no_removed <- length(removed)
    no_staying <- length(actualprey)
    
    # individuals in actualprey reproduce and produce no_removed offspring, to keep population size constant
    if (haploid_asexual == T) { ### HAPLOID ASEXUAL
      alleles[,t] <- alleles[,t-1]
      
      for (i in 1:no_removed) {
        k <- sample(1:no_staying, 1)
        ii <- actualprey[k]
        jj <- removed[i]
        age[jj] <- 0
        alleles[jj,t] <- alleles[ii,t-1] # inheritance
      }
    }
    else { ### HAPLOID SEXUAL
      alleles[,t] <- alleles[,t-1]

      dalleles <- rep(0, no_removed)
      for (i in 1:no_removed) {
        mates <- sample(actualprey, 2, replace=F)
        if (runif(1) < 0.5) {
          dalleles[i] <- alleles[mates[1],t-1]
        }
        else {
          dalleles[i] <- alleles[mates[2],t-1]
        }
        jj <- removed[i]
        age[jj] <- 0
        alleles[jj,t] <- dalleles[i] # inheritance
        escape_prob_temp[jj] <- init_escape_prob
      }
    }
    
    if (immediate_update == T) { # within timestep update
      escape_prob[,t] <- escape_prob[,t-1] 
    }
    else { # between timestep update
      escape_prob[,t] <- escape_prob_temp 
    }
    escape_prob[removed,t] <- rep(init_escape_prob, no_removed) # new individuals have no escape experience
    
    # calculate number of unique alleles
    no_diff_alleles[s,t] <- length(unique(alleles[,t]))

  } # loop over t
  #allelesend[,s] <- alleles[,maxtime]
  
} # loop over s

### without learning ###

#allelesbeg2 <- array(0, c(no_individuals, no_simulations))
#allelesend2 <- array(0, c(no_individuals, no_simulations))
no_diff_alleles2 <- array(0, c(no_simulations, maxtime))
#totalremoved2 <- rep(0, no_simulations)
#totalattacks2 <- rep(0, no_simulations)
age2 <- rep(0, no_individuals)
#totalremoved2x <- rep(0, no_simulations)  
age_at_death2 <- NULL
age_dist2 <- array(0, c(age_limit, no_simulations))

for (s in 1:no_simulations) {
  
  age2 <- sample(0:(lifespan-1), no_individuals, replace=T)
  alleles2 <- array(0, c(no_individuals, maxtime))
  alleles2[,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  #allelesbeg2[,s] <- alleles2[,1]
  no_diff_alleles2[s,1] <- length(unique(alleles2[,1]))
  escape_prob2 <- array(init_escape_prob, c(no_individuals, maxtime))
  
  for (t in 2:maxtime) {
    # saving age distribution
    if (t == age_time_sample) {
      for (i in 1:(age_limit-1)) {
        age_dist2[i,s] <- sum(age2 == i-1)
      }
      age_dist2[age_limit,s] <- sum(age2 >= (age_limit-1))
      show(sum(age_dist2))
    }
    
    # aging
    age2 <- age2 + 1

    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed2 <- NULL
    
    actualprey2 <- 1:no_individuals
    repeat {
      #totalattacks2[s] <- totalattacks2[s] + 1
      attack <- sample(actualprey2, 1) # random choice of prey
      if (runif(1) >= escape_prob2[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed2 <- c(attack, removed2)
        actualprey2 <- setdiff(actualprey2, attack)
        age_at_death2 <- c(age_at_death2, age2[attack])
      } 
      if (successful_attacks >= no_remove) { break }
    }
    #show(successful_attacks)
    #show(length(removed))
    #totalremoved2[s] <- totalremoved2[s] + length(removed)
    
    # remove surviving old individuals
    if (aging) {
      apC2 <- length(actualprey2)
      if(!random_aging) {
        for (i in apC2:-1:1) {
          ii <- actualprey2[i]
          if (age2[ii] > lifespan) {
            age_at_death2 <- c(age_at_death2, age2[ii])
            removed2 <- c(ii, removed2)
            actualprey2 <- setdiff(actualprey2, ii)
          }
        }
      }
      else { # random_aging == T
        for (i in apC2:-1:1) {
          ii <- actualprey2[i]
          if (runif(1) < 1/lifespan) {
            age_at_death2 <- c(age_at_death2, age2[ii])
            removed2 <- c(ii, removed2)
            actualprey2 <- setdiff(actualprey2, ii)
          }
        }
      }
      #totalremoved2x[s] <- totalremoved2x[s] + length(removed2)
      #staying2 <- setdiff(staying, removed2) # finally staying
      #removed <- c(removed, removed2) # total removed
    }
    no_removed <- length(removed2)
    no_staying <- length(actualprey2)
    
    # individuals in actualprey2 reproduce and produce no_remove offspring, to keep population size constant
    if (haploid_asexual == T) { ### HAPLOID ASEXUAL
      alleles2[,t] <- alleles2[,t-1]
      
      for (i in 1:no_removed) {
        k <- sample(1:no_staying, 1)
        ii <- actualprey2[k]
        jj <- removed2[i]
        age2[jj] <- 0
        alleles2[jj,t] <- alleles2[ii,t-1] # inheritance
      }
    }
    else { ### HAPLOID SEXUAL
      alleles2[,t] <- alleles2[,t-1]
      
      dalleles <- rep(0, no_removed)
      for (i in 1:no_removed) {
        mates <- sample(actualprey2, 2, replace=F)
        if (runif(1) < 0.5) {
          dalleles[i] <- alleles2[mates[1],t-1]
        }
        else {
          dalleles[i] <- alleles2[mates[2],t-1]
        }
        jj <- removed2[i]
        age2[jj] <- 0
        alleles2[jj,t] <- dalleles[i] # inheritance
      }
    }
  
    # calculate number of unique alleles
    no_diff_alleles2[s,t] <- length(unique(alleles2[,t]))

  } # loop over t
  #allelesend2[,s] <- alleles2[,maxtime]

} # loop over s

library("RColorBrewer")
col <- brewer.pal(n = 3, name = "Paired")

ccex <- 1.5
ccex2 <- 1.5

#par(mfrow=c(1,1))
xmax <- maxtime
ymax <- init_no_alleles
plot(1:maxtime, no_diff_alleles[1,], type = 'l', lwd = 3, col = col[1], xlab = 'Time', ylab = 'Number of different alleles', cex.axis = ccex, cex.lab = ccex, xlim = c(1, maxtime), ylim = c(0, ymax))
lines(1:maxtime, no_diff_alleles2[1,], type = 'l', lwd = 3, col = col[2])
legend('topright', c('with learning', 'without learning'), lty = 1, pch = -1, bty = 'o', lwd = 3, cex = ccex2, col = col)
for (s in 2:no_simulations) {
  lines(1:maxtime, no_diff_alleles[s,], type = 'l', lwd = 3, col = col[1])
  lines(1:maxtime, no_diff_alleles2[s,], type = 'l', lwd = 3, col = col[2])
}

# png("hist.png", width=800, height=600)
# par(mfrow=c(2,2))
# ss <- 2
# hist(allelesbeg[,ss], main="With learning: t=0", xlab="Allele number", breaks=0:100, xlim=c(1,100))
# hist(allelesend[,ss], main="With learning: t=100", xlab="Allele number", breaks=0:100, xlim=c(1,100))
# hist(allelesbeg2[,ss], main="Without learning: t=0", xlab="Allele number", breaks=0:100, xlim=c(1,100))
# hist(allelesend2[,ss], main="Without learning: t=100", xlab="Allele number", breaks=0:100, xlim=c(1,100))
# par(mfrow=c(1,1))
# dev.off()

# show(totalremoved/(maxtime-1))
# show(totalremoved2/(maxtime-1))
# show(totalremovedx/(maxtime-1))
# show(totalremoved2x/(maxtime-1))
# show(totalattacks/(maxtime-1))
# show(totalattacks2/(maxtime-1))

# png("hist.png", width=800, height=600)
#par(mfrow=c(2,2))

# means
bilan <- array(0,c(2,age_limit))
bilan[1,] <- rowMeans(age_dist)
bilan[2,] <- rowMeans(age_dist2)
rownames(bilan) <- c("with learning", "without learning")
colnames(bilan) <- 0:(age_limit-1)
bilan <- as.matrix(bilan)

#Plot boundaries
lim <- 1.2*max(bilan)

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# sds
library(matrixStats)
stdev <- array(0,c(2,age_limit))
stdev[1,] <- rowSds(age_dist)
stdev[2,] <- rowSds(age_dist2)
rownames(stdev) <- c("with learning", "without learning")
colnames(stdev) <- 0:(age_limit-1)
stdev <- as.matrix(stdev)

#I am ready to add the error bar on the plot using my "error bar" function !
ze_barplot <- barplot(bilan, beside=T, legend.text=F, col=c("blue", "skyblue"), ylim=c(0,lim), xlab="Age", ylab="Number of individuals", cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
error.bar(ze_barplot, bilan, stdev)
legend("topright", 
       legend = c("with learning", "without learning"), 
       fill = c("blue", "skyblue"), cex = 1.5)
#par(mfrow=c(1,1))
#dev.off()

show(mean(age_at_death))
show(sd(age_at_death))
show(mean(age_at_death2))
show(sd(age_at_death2))
