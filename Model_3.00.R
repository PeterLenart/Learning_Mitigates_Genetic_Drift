# working directory
setwd("c:/WORK/Rozpracovano/Genetic drift and predation/")

############
# parameters
############
no_individuals <- 1000 # population size
init_no_alleles <- 100 # initial number of unique alleles
maxtime <- 100 # simulation time
no_remove <- 500 # number of individuals (prey) consumed each timestep
init_escape_prob <- 0.5 # initial probability of escaping predator attack
max_escape_prob <- 0.99 # maximum probability of escaping predator attack
rate_escape_prob <- 0.75 # proportional increase in escape probability (for learning; 0 for no learning)
no_simulations <- 5 # number of simulation runs
haploid_asexual <- T # haploid asexual (T) or haploid sexual (F) prey
aging <- F # finite lifespan set on (T) or off (F)
lifespan <- 10 # maximum age if aging = T
immediate_update <- F # original update once per time step (F) or after every unsuccessful attack
############

### with learning ###

allelesbeg <- array(0, c(no_individuals, no_simulations))
allelesend <- array(0, c(no_individuals, no_simulations))
no_diff_alleles <- array(0, c(no_simulations, maxtime))
totalremoved <- rep(0, no_simulations)
totalattacks <- rep(0, no_simulations)
age <- rep(0, no_individuals)
totalremovedx <- rep(0, no_simulations)
age_at_death <- NULL

for (s in 1:no_simulations) {
  
  escape_prob <- array(init_escape_prob, c(no_individuals, maxtime))
  alleles <- array(0, c(no_individuals, maxtime))
  alleles[,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  allelesbeg[,s] <- alleles[,1]
  no_diff_alleles[s,1] <- length(unique(alleles[,1]))
  age <- sample(0:lifespan, no_individuals, replace=T)
  
  for (t in 2:maxtime) {
    # aging
    if (t > 2) age <- age + 1
    
    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    actualprey <- 1:no_individuals

    repeat {
      attack <- sample(actualprey, 1) # random choice of prey
      totalattacks[s] <- totalattacks[s] + 1
      if (runif(1) >= escape_prob[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack, removed)
        actualprey <- setdiff(actualprey, attack)
        age_at_death <- c(age_at_death, age[attack])
      } else { # unsuccessful predator attack: increase in prey escape probability with learning
        if (immediate_update == F) { # between timestep update
          escape_prob[attack,t] <- min(max_escape_prob,(1+rate_escape_prob)*escape_prob[attack,t-1])
        }
        else { # within timestep update
          escape_prob[attack,t-1] <- min(max_escape_prob,(1+rate_escape_prob)*escape_prob[attack,t-1]) 
        }
      }
      if (successful_attacks >= no_remove) { break }
    }
    #show(successful_attacks)
    #show(length(removed))
    totalremoved[s] <- totalremoved[s] + length(removed)
    
    # determining indices of non-removed individuals
    staying <- setdiff(1:no_individuals, removed)
    #show(length(staying))
    
    # remove surviving old individuals
    staying2 <- staying
    if (aging) {
      removed2 <- NULL
      for (i in staying) {
        if (age[i] > lifespan) {
          removed2 <- c(i, removed2)
          age_at_death <- c(age_at_death, age[i])
        }
      }
      totalremovedx[s] <- totalremovedx[s] + length(removed2)
      staying2 <- setdiff(staying, removed2) # finally staying
      removed <- c(removed, removed2) # total removed
    }
    no_removed <- length(removed)
    
    # staying2 individuals reproduce and produce no_removed offspring, to keep population size constant
    if (haploid_asexual == T) { ### HAPLOID ASEXUAL
      reproducing <- sample(staying2, no_removed, replace=T)
      alleles[,t] <- alleles[,t-1]
      alleles[removed,t] <- alleles[reproducing,t]
      age[removed] <- rep(0, no_removed)
    }
    else { ### HAPLOID SEXUAL
      dalleles <- vector(mode="double", length=no_removed)
      for (i in 1:no_removed) {
        mates <- sample(staying2, 2, replace=F)
        if (runif(1) < 0.5) {
          dalleles[i] <- alleles[mates[1],t-1]
        }
        else {
          dalleles[i] <- alleles[mates[2],t-1]
        }        
      }
      alleles[,t] <- alleles[,t-1]
      alleles[removed,t] <- dalleles
      age[removed] <- rep(0, no_removed)
    }
    if (immediate_update == T) { # within timestep update
      escape_prob[,t] <- escape_prob[,t-1] 
    }
    escape_prob[removed,t] <- rep(init_escape_prob, no_removed) # new individuals have no escape experience

    # calculate number of unique alleles
    no_diff_alleles[s,t] <- length(unique(alleles[,t]))

  } # loop over t
  allelesend[,s] <- alleles[,maxtime]
  
} # loop over s

### without learning ###

allelesbeg2 <- array(0, c(no_individuals, no_simulations))
allelesend2 <- array(0, c(no_individuals, no_simulations))
no_diff_alleles2 <- array(0, c(no_simulations, maxtime))
totalremoved2 <- rep(0, no_simulations)
totalattacks2 <- rep(0, no_simulations)
age <- rep(0, no_individuals)
totalremoved2x <- rep(0, no_simulations)  
age_at_death2 <- NULL

for (s in 1:no_simulations) {
  
  escape_prob2 <- array(init_escape_prob, c(no_individuals, maxtime))
  alleles <- array(0, c(no_individuals, maxtime))
  alleles[,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  allelesbeg2[,s] <- alleles[,1]
  no_diff_alleles2[s,1] <- length(unique(alleles[,1]))
  age <- sample(0:lifespan, no_individuals, replace=T)
  
  for (t in 2:maxtime) {
    # aging
    if (t > 2) age <- age + 1

    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    actualprey <- 1:no_individuals

    repeat {
      totalattacks2[s] <- totalattacks2[s] + 1
      attack <- sample(actualprey, 1) # random choice of prey
      if (runif(1) >= escape_prob2[attack,t-1]) { # successful predator attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack, removed)
        actualprey <- setdiff(actualprey, attack)
        age_at_death2 <- c(age_at_death2, age[attack])
      } 
      if (successful_attacks >= no_remove) { break }
    }
    #show(successful_attacks)
    #show(length(removed))
    totalremoved2[s] <- totalremoved2[s] + length(removed)
    
    # determining indices of non-removed individuals
    staying <- setdiff(1:no_individuals, removed)
    #show(length(staying))
    
    # remove surviving old individuals
    staying2 <- staying
    if (aging) {
      removed2 <- NULL
      for (i in staying) {
        if (age[i] > lifespan) {
          removed2 <- c(i, removed2)
          age_at_death2 <- c(age_at_death2, age[i])
        }
      }
      totalremoved2x[s] <- totalremoved2x[s] + length(removed2)
      staying2 <- setdiff(staying, removed2) # finally staying
      removed <- c(removed, removed2) # total removed
    }
    no_removed <- length(removed)
    
    # staying individuals reproduce and produce no_remove offspring, to keep population size constant
    if (haploid_asexual == T) { ### HAPLOID ASEXUAL
      reproducing <- sample(staying2, no_removed, replace=T)
      alleles[,t] <- alleles[,t-1]
      alleles[removed,t] <- alleles[reproducing,t]
      age[removed] <- rep(0, no_removed)
    }
    else { ### HAPLOID SEXUAL
      dalleles <- vector(mode="double", length=no_removed)
      for (i in 1:no_removed) {
        mates <- sample(staying2, 2, replace=F)
        if (runif(1) < 0.5) {
          dalleles[i] <- alleles[mates[1],t-1]
        }
        else {
          dalleles[i] <- alleles[mates[2],t-1]
        }        
      }
      alleles[,t] <- alleles[,t-1]
      alleles[removed,t] <- dalleles
      age[removed] <- rep(0, no_removed)
    }
    
    # calculate number of unique alleles
    no_diff_alleles2[s,t] <- length(unique(alleles[,t]))

  } # loop over t
  allelesend2[,s] <- alleles[,maxtime]

} # loop over s

library("RColorBrewer")
col <- brewer.pal(n = 3, name = "Paired")

ccex <- 1.5
ccex2 <- 1.5

par(mfrow=c(1,1))
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

show(mean(age_at_death))
show(sd(age_at_death))
show(mean(age_at_death2))
show(sd(age_at_death2))
