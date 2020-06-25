
no_individuals <- 1000
init_no_alleles <- 100 # initial number of unique alleles
maxtime <- 100
no_remove <- 800 # number of individuals consumed each timestep
init_escape_prob <- 0.5 # initial probability of escaping attack
no_simulations <- 10
maxescape <- 0.95

### with learning

rate_escape_prob <- 0.75 # proportional increase in escape probability

no_diff_alleles <- array(0, c(no_simulations, maxtime))

for (s in 1:no_simulations) {
  
  escape <- array(init_escape_prob, c(no_individuals, maxtime))
  
  alleles <- array(0, c(no_individuals, maxtime))
  alleles[,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  
  no_diff_alleles[s,1] <- length(unique(alleles[,1]))
  
  for (t in 2:maxtime) {
    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    
    repeat {
      attack <- sample(1:no_individuals, 1) # random choice of prey
      if (runif(1) >= escape[attack,t-1]) { # successful attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack,removed)
      } else { # unsuccessful attack: increase in escape probability with learning
        escape[attack,t] <- min(maxescape,(1+rate_escape_prob)*escape[attack,t-1])
      }
      if (successful_attacks >= no_remove) { break }
    }

    # determining indices of not removed individuals
    b <- 1:no_individuals
    staying <- b [! b %in% removed]
    
    # no_remove of the remaining individuals reproduce, to keep population size constant
    reproducing <- sample(staying, no_remove, replace=T)
    alleles[,t] <- alleles[,t-1]
    alleles[removed,t] <- alleles[reproducing,t]
    escape[removed,t] <- rep(init_escape_prob, no_remove) # new individuals have no escape experience
      
    # calculate number of unique alleles
    no_diff_alleles[s,t] <- length(unique(alleles[,t]))

  } # loop over t

} # loop over s

### without learning (same as above, but with rate_escape_prob = 0 [meaning no increase])

rate_escape_prob <- 0

no_diff_alleles2 <- array(0, c(no_simulations, maxtime))
  
for (s in 1:no_simulations) {
  
  escape <- array(init_escape_prob, c(no_individuals, maxtime))
  
  alleles <- array(0, c(no_individuals, maxtime))
  alleles[,1] <- sample(1:init_no_alleles, no_individuals, replace=T)
  
  no_diff_alleles2[s,1] <- length(unique(alleles[,1]))
  
  for (t in 2:maxtime) {
    # remove no_remove individuals: loop over attacks
    successful_attacks <- 0
    removed <- NULL
    
    repeat {
      attack <- sample(1:no_individuals, 1)
      if (runif(1) >= escape[attack,t-1]) { # successful attack
        successful_attacks <- successful_attacks + 1
        removed <- c(attack,removed)
      } else { # unsuccessful attack: increase in escape probability with learning
        escape[attack,t] <- min(maxescape,(1+rate_escape_prob)*escape[attack,t-1])
      }
      if (successful_attacks >= no_remove) { break }
    }
    
    b <- 1:no_individuals
    staying <- b [! b %in% removed]
    
    # no_remove of the remaining individuals reproduce
    reproducing <- sample(staying, no_remove, replace=T)
    alleles[,t] <- alleles[,t-1]
    alleles[removed,t] <- alleles[reproducing,t]
    
    # calculate number of unique alleles
    no_diff_alleles2[s,t] <- length(unique(alleles[,t]))

  } # loop over t

} # loop over s

library("RColorBrewer")
col = brewer.pal(n = 3, name = "Paired")

ccex <- 1.5
ccex2 <- 1.5

xmax <- maxtime
ymax <- init_no_alleles
plot(1:maxtime, no_diff_alleles[1,], type = 'l', lwd = 3, col = col[1], xlab = 'Time', ylab = 'Number of different alleles', cex.axis = ccex, cex.lab = ccex, xlim = c(1, xmax), ylim = c(0, ymax))
lines(1:maxtime, no_diff_alleles2[1,], type = 'l', lwd = 3, col = col[2])
legend('bottomleft', c('with learning', 'without learning'), lty = 1, pch = -1, bty = 'o', lwd = 3, cex = ccex2, col = col)
for (s in 2:no_simulations) {
  lines(1:maxtime, no_diff_alleles[s,], type = 'l', lwd = 3, col = col[1])
  lines(1:maxtime, no_diff_alleles2[s,], type = 'l', lwd = 3, col = col[2])
}
