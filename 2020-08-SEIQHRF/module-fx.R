##
## SEIR Model: Adding an Exposed State to an SIR
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Venkata R. Duvvuri
## Date: August 2018
##

# Replacement discord_edgelist function -----------------------------------

# Add infective_status parameters to the original discord_edgelist function
# So more than one status can be infective 

discord_edgelist <- function (dat, at, inf_status, network = 1) {
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
        isInf <- matrix(stat %in% inf_status, ncol = 2)
        isSus <- matrix(stat %in% "s", ncol = 2)
        SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
        ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
        pairs <- rbind(SIpairs, ISpairs[, 2:1])
        if (nrow(pairs) > 0) {
            sus <- pairs[, 1]
            inf <- pairs[, 2]
            inf_status <- status[inf]
            del <- data.frame(at, sus, inf, inf_status)
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

  # Uncomment this to function environment interactively
  # browser()

  # Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  # Parameters 
  infections <- get_param(dat, "infections")
  n_infections <- nrow(infections)
  infective_status <- infections$from 
  infective_ids <- which(active == 1 & status %in% infective_status)

  # Retroactively save infection time as 1 for ids initialized in
  # infective status.
  if (at == 2) {
      infTime <- rep(NA, length(active))
      infTime[infective_ids] <- 1
      dat <- set_attr(dat, "inf_time", infTime)
  } else {
      infTime <- get_attr(dat, "inf_time")
  }

  del <- discord_edgelist(dat, at, infective_status)

  for (i in 1:n_infections) {

      del_i <- del[del$status == infective_status[i], ]

      # Stochastic transmission process
      transmit <- rbinom(nrow(del_i), 1, infections$final.prob[i])

      # Keep rows where transmission occurred
      del_i <- del_i[which(transmit == 1), ]

      # Look up new ids if any transmissions occurred
      idsNewInf <- unique(del_i$sus)

      # Set new attributes for those newly infected
      status[idsNewInf] <- "e"
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "status", status)
      dat <- set_attr(dat, "infTime", infTime)

      ## Save summary statistic for S->E flow
      # TODO decide about this
      # dat <- set_epi(dat, "se.flow", at, nInf)

  }

  return(dat)

}

# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  # Uncomment this to function environment interactively
  # browser()

  # Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  # Progressions
  progressions <- get_param(dat, "progressions")
  n_progressions <- length(progressions)
  from <- progressions$from
  to <- progressions$to
  rate <- progressions$rate

  for (p in 1:n_progressions) {

    from_p <- from[p]
    to_p <- to[p]
    rate_p <- rate[p]

    # Sample progression ids
    n_eligible <- sum(active == 1 & status == from_p)
    ids <- which(rbinom(n_eligible, 1, rate_p) == 1)

    # Save flow statistic
    dat <- set_epi(dat, paste0(from_p, to_p, ".flow"), at, length(ids))

    # Update status
    status[ids] <- to_p
  
  }

  # Save updated status
  dat <- set_attr(dat, "status", status)

  # Save num statistic
  for (s in unique(c(from, to))) {

    dat <- set_epi(dat, paste0(s, ".num"), at, sum(active == 1 & status == s))
  
  }

  return(dat)

}

