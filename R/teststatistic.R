###############################
# Compute Test Statistic      #
###############################

#===========================================================================

# input:  a data set of the function prepare
# output: the test statistic and p-value

#===========================================================================


# KM  ##########################################################################

#-------------------------------------------------------------------------------
# input:  a data set of the function prepare
# output: - Kaplan-Meier estimates (left- and right-continuous version)
#         - Maximum of the event times
#         - Number at risk (left- and right-continuous version)
#         - Number of events (left- and right-continuous version)
#         - empirical distribution function (left- and right-continuous version)

#           this is made for the whole data set as well as for the subsets
#           returns the values at all event times
#-------------------------------------------------------------------------------


KM <- function(data) {

  # for the whole data set ----------------------------------------------------

  n <- dim(data[[1]])[1]

  surv_all <- survival::survfit(Surv(data[[1]][, 1], data[[1]][, 2]) ~ 1)

  # Kaplan-Meier estimates
  km_help_left  <- stepfun(surv_all$time, c(1, surv_all$surv), right=TRUE)
  km_help_right <- stepfun(surv_all$time, c(1, surv_all$surv), right=FALSE)
  km            <- km_help_right(surv_all$time)
  km_minus      <- km_help_left(surv_all$time )

  # Maximum of all event times
  max <- NULL

  # Number at risk
  Y_help_left  <- stepfun(surv_all$time, c(surv_all$n.risk,surv_all$n.censor[length(surv_all$n.censor)]), right=TRUE)
  Y_help_right <- stepfun(surv_all$time, c(surv_all$n.risk,surv_all$n.censor[length(surv_all$n.censor)]), right=FALSE)
  Y            <- Y_help_right(surv_all$time)
  Y_minus      <- Y_help_left(surv_all$time)

  # Number of events
  N_help_left  <- stepfun(surv_all$time, cumsum(c(0,surv_all$n.event)), right=TRUE)
  N_help_right <- stepfun(surv_all$time, cumsum(c(0,surv_all$n.event)), right=FALSE)
  N            <- N_help_right(surv_all$time)
  N_minus      <- N_help_left(surv_all$time)


  # Kaplan-Meier estimator of distribution function
  F_          <- 1 - 0.5*(km+km_minus)
  F_minus     <- c(0, F_[1:(length(F_) - 1)])

  # for the subsets ------------------------------------------------------------

  help <- length(levels(data[[1]]$cell))

  for (i in 1:help) {

    n    <- c(n, dim(data[[i+1]])[1])

    surv <- survival::survfit(Surv(data[[i+1]][, 1], data[[i+1]][, 2]) ~ 1)

    # Kaplan-Meier estimates
    km_help_left  <- stepfun(surv$time, c(1, surv$surv), right=TRUE)
    km_help_right <- stepfun(surv$time, c(1, surv$surv), right=FALSE)
    km            <- cbind(km, km_help_right(surv_all$time))
    km_minus      <- cbind(km_minus, km_help_left(surv_all$time))

    # Maximum of the event times
    max           <- c(max, max(data[[i+1]][1][data[[i+1]][2] == 1]))

    # Number at risk
    Y_help_left  <- stepfun(surv$time, c(surv$n.risk, surv$n.censor[length(surv$n.censor)]), right=TRUE)
    Y_help_right <- stepfun(surv$time, c(surv$n.risk, surv$n.censor[length(surv$n.censor)]), right=FALSE)
    Y            <- cbind(Y, Y_help_right(surv_all$time))
    Y_minus      <- cbind(Y_minus, Y_help_left(surv_all$time ))

    # Number of events
    N_help_left  <- stepfun(surv$time, cumsum(c(0,surv$n.event)), right=TRUE)
    N_help_right <- stepfun(surv$time, cumsum(c(0,surv$n.event)), right=FALSE)
    N            <- cbind(N, N_help_right(surv_all$time))
    N_minus      <- cbind(N_minus, N_help_left(surv_all$time))

    # Kaplan-Meier estimator of distribution function
    F_help_right <- 1 - 0.5*(km_help_right(surv_all$time) + km_help_left(surv_all$time))
    F_help_left  <- c(0, F_help_right[1:(length(F_help_right) - 1)])
    F_           <- cbind(F_, F_help_right)
    F_minus      <- cbind(F_minus, F_help_left)
  }


  result <- c(list(n), list(surv_all$time), list(max), list(Y), list(Y_minus), list(N),
              list(N_minus) ,list(km), list(km_minus), list(F_), list(F_minus))

  return (result)
}




# S_H_hat ##########################################################################

#-------------------------------------------------------------------------------
# input:  an object of function KM
# output: values of function S_H_hat at all event times
#-------------------------------------------------------------------------------

S_H_hat <- function(km) {
  return ((0.5) * (1 / km[[1]][1]) * (km[[4]][, 1] + km[[5]][, 1]))
}



# integral ##########################################################################

#-------------------------------------------------------------------------------
# input:  - an object of function KM
#         - an object of function S_H_hat
# output: values of integral for all subsets and all event times
#-------------------------------------------------------------------------------


integral <- function(km,S_H_hat){

  T_     <- round(min(km[[3]]),7)
  ind    <- which(round(km[[2]],7) == T_)
  help   <- (S_H_hat * (km[[10]] - km[[11]]))
  result <- apply(X = matrix(help[ind, ], 1), MARGIN = 2, FUN = sum)

  for (i in (ind - 1):1) {
    help1  <- apply(X = help[i:ind, ], MARGIN = 2, FUN = sum)
    result <- rbind(help1, result)

  }

  return(result)

}


# h_hat ##########################################################################

#-------------------------------------------------------------------------------
# input:  - an object of function KM
#         - an object of function S_H_hat
#         - an object of function integral
# output: values of h_hat for all subsets and all event times
#-------------------------------------------------------------------------------



h_hat <- function (km, S_H_hat, integral) {
  T_     <- round(min(km[[3]]),7)
  ind    <- which(round(km[[2]],7) == T_)
  result <- (km[[9]][1:ind, ])^2 * (S_H_hat[1:ind] - ((1 / km[[8]][1:ind, ]) * integral))^2 * (1 / ((1 / km[[1]]) * km[[5]][1:ind, ]))
  return (result)
}



# sigma ##########################################################################

#-------------------------------------------------------------------------------
# input:  - an object of function KM
#         - an object of function h_hat
# output: values of sigma for all subsets and all event times
#-------------------------------------------------------------------------------


sigma <- function(km,h_hat) {

  T_     <- round(min(km[[3]]),7)
  ind    <- which(round(km[[2]],7) == T_)
  help   <- h_hat * (1 - (((km[[6]][1:ind, ] - km[[7]][1:ind, ] - 1)) / (km[[5]][1:ind, ] - 1))) * (1/km[[5]][1:ind, ]) *
            (km[[6]][1:ind, ] - km[[7]][1:ind, ])

  result   <- apply(X = help[1:ind, ], MARGIN = 2, FUN = sum)



  return(result)
}



# statistic ##########################################################################

#-------------------------------------------------------------------------------
# input:  - an object of function KM
#         - an object of function sigma
#         - an object of function integral
#         - an object of function C
# output: value of the test statistic
#-------------------------------------------------------------------------------



statistic <- function (km, sigma, integral, C) {

  help1   <- length(km[[3]])+1
  help    <- C %*% integral[1, ][2:help1]
  help2   <- (km[[1]][1] / km[[1]][2:help1]) * sigma[2:help1]
  V_hat   <- diag(help2)
  stat    <- km[[1]][1] * t(help) %*% solve(C %*% V_hat %*% t(C)) %*% help
  k       <- qr(C)$rank
  p_value <-  1- pchisq(stat, k)

  return (c(stat, p_value))

}








