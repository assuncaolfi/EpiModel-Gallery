##
## SEIR Model: Adding an Exposed State to an SIR
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Venkata R. Duvvuri
## Date: August 2018
##


# Replacement discord_edgelist function -----------------------------------

# Add infective_status parameters to the original discord_edgelist function
# In our case, some status are infective

discord_edgelist <- function (dat, at, infectiveStatus, network = 1) {
    status <- get_attr(dat, "status")
    active <- get_attr(dat, "active")
    tergmLite <- get_control(dat, "tergmLite")
    if (tergmLite == TRUE) {
        el <- dat$el[[network]]
    }
    else {
        el <- get.dyads.active(dat$nw[[network]], at = at)
    }
    del <- NULL
    if (nrow(el) > 0) {
        el <- el[sample(1:nrow(el)), , drop = FALSE]
        stat <- matrix(status[el], ncol = 2)
        isInf <- matrix(stat %in% infectiveStatus, ncol = 2)
        isSus <- matrix(stat %in% "s", ncol = 2)
        SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
        ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
        pairs <- rbind(SIpairs, ISpairs[, 2:1])
        if (nrow(pairs) > 0) {
            sus <- pairs[, 1]
            inf <- pairs[, 2]
            from <- status[inf]
            del <- data.frame(at, sus, inf, from)
            keep <- rowSums(matrix(c(active[del$sus], active[del$inf]),
                ncol = 2)) == 2
            del <- del[keep, ]
            if (nrow(del) < 1) {
                del <- NULL
            }
        }
    }
    return(del)
}

# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  if (at == 2) {
    infTime <- rep(NA, length(active))
    infTime[which(status == "i")] <- 1
    dat <- set_attr(dat, "infTime", infTime)
  } else {
    infTime <- get_attr(dat, "infTime")
  }

  ## Parameters ##
  inf.pars <- get_param(dat, "inf.pars")
  prog.pars <- get_param(dat, "prog.pars")
  sum.pars <- get_param(dat, "sum.pars")

  ## Find infected nodes ##
  infectedStatus <- sum.pars$infected.status
  idsInf <- which(active == 1 & status %in% infectedStatus)
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize default incidence at 0 ##
  nInf <- 0

  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {

    ## Look up discordant edgelist ##
    infectiveStatus <- sum.pars$infective.status
    del <- discord_edgelist(dat, at, infectiveStatus)

    ## If any discordant pairs, proceed ##
    if (!is.null(del)) {

      # Set parameters on discordant edgelist data frame
      del <- merge(del, inf.pars, by = "from")

      # Stochastic transmission process
      transmit <- rbinom(nrow(del), 1, del$final.prob)

      # Keep rows where transmission occurred
      del <- del[which(transmit == 1), ]

      # Look up new ids if any transmissions occurred
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)

      # Set new attributes for those newly infected
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  ## Save summary statistic for S->E flow
  dat <- set_epi(dat, "se.flow", at, nInf)

  return(dat)
}

# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Uncomment this to function environment interactively
  # browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  previousStatus <- status

  ## Parameters ##
  progPars <- get_param(dat, "prog.pars")
  sum.pars <- get_param(dat, "sum.pars")

  ## X to Y progression process ##
#   status <- sample(c(prog.pars$from, "s"), 20, replace = TRUE)
  progRates <- sum.pars$prog.rates
  probProg <- progRates[status]
  probProg[is.na(probProg)] <- 0
  isProg <- active & rbinom(length(status), 1, probProg)
  noProg <- all(isProg == FALSE)

  if (noProg) {
    status <- status
  } else {
    nextStates <- sum.pars$next.states
    nextProbs <- sum.pars$next.probs
    statusProg <- status[isProg == 1]
    status[isProg == 1] <- mapply(
      sample,
      x = nextStates[statusProg],
      size = 1,
      prob = nextProbs[statusProg]
    )
  }

  ## Write out updated status attribute ##
  dat <- set_attr(dat, "status", status)

  ## Save flows ##
  uniqueFlowNames <- sum.pars$flow.names
  vecFlows <- paste0(previousStatus, status, ".flow")[isProg]
  progFlows <- table(vecFlows)
  progFlowNames <- names(progFlows)
  for (flowName in uniqueFlowNames) {
    value <- if (flowName %in% progFlowNames) progFlows[flowName] else 0
    dat <- set_epi(dat, flowName, at, value)
  } 

  ## Save nums ##
  uniqueNumNames <- sum.pars$num.names
  progNums <- table(status[active == 1])
  progNumNames <- paste0(names(progNums), ".num")
  for (numName in uniqueNumNames) {
    value <- if (numName %in% progNumNames) progNums[numName] else 0
    dat <- set_epi(dat, numName, at, value)
  }

  return(dat)
}

