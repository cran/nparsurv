#' Nonparametric Tests for Main Effects, Simple Effects and Interaction Effect in a Two-Factorial Design with Censored Data
#'
#'The \code{nparsurv_test} function calculates the test statistics and the p-values as described in 'Nonparametric Methods for Factorial Designs with Censored Data'
#'by Akritas and Brunner.
#'
#'
#' @param data A \code{\link{data.frame}} of the form (time, status, factorA, factorB)
#'    \describe{
#'    \item{time:}{time of event or censoring, numeric}
#'     \item{status:}{indicator for censoring, 1=event, 0=censored, integer}
#'    \item{factorA:}{first factor variable, factor}
#'      \item{factorB:}{second factor variable, factor}
#'    }
#'
#'
#' Missing values must be saved as \code{NA}.
#'
#' @details The package provides tests for a survival setting with two influencing variables, that are factors with at
#'  least two levels each. Details are shown in 'Nonparametric Methods for Factorial Designs with Censored Data' by Akritas
#'   and Brunner.
#'   The \code{nparsurv_test} function returns the values of the five test statistics: the tests for main effects, simple effects
#'   and the interaction effect. Additionally, based on the asymptotic chi-square distribution of the test statistic under the nullhypothesis, p-values are computed.
#'
#'
#' @return A \code{nparsurv_test} object containing the following components:
#' \item{maineffectA / maineffectB}{The test statistic and p-value for the nullhypotheses 'no main effect of factor A' and 'no main effect of factor B' respectively.}
#'  \item{simpleeffectA / simpleeffectB}{The test statistic and p-value for the null hypotheses 'no simple effect of factor A' and 'no simple effect of factor B' respectively.}
#'  \item{interactioneffect}{The test statistic and p-value for the null hypothesis 'no interaction effect between factor A and factor B'.}
#'
#' @examples
#' data_ovarian<-data.frame(survival::ovarian$futime,
#'                        survival::ovarian$fustat,
#'                        as.factor(survival::ovarian$resid.ds),
#'                        as.factor(survival::ovarian$rx))
#' nparsurv_test(data_ovarian)
#'
#' data_GBSG2<-data.frame(TH.data::GBSG2$time,
#'                        TH.data::GBSG2$cens,
#'                        TH.data::GBSG2$tgrade,
#'                        TH.data::GBSG2$horTh)
#' nparsurv_test(data_GBSG2)
#'
#' @references Michael G. Akritas, Edgar Brunner(1997). Nonparametric Methods for Factorial Designs with Censored Data. Journal
#'            of the American Statistical Association.
#' @import survival
#' @import TH.data
#' @importFrom stats na.omit pchisq stepfun
#'
#' @export


nparsurv_test <- function (data) {


  levels3    <- levels(data[, 3])
  a          <- length(levels3)
  levels4    <- levels(data[, 4])
  b          <- length(levels4)
  help1      <- rep(levels3, each = b)
  help2      <- rep(levels4, a)
  cells      <- matrix(c(help1,help2), , 2)
  C          <- contrast(a, b, c(1, 1, 1, 1, 1))
  data       <- prepare(data)
  km         <- KM(data)
  S_H_hat    <- S_H_hat(km)
  integral   <- integral(km, S_H_hat)
  h_hat      <- h_hat(km, S_H_hat, integral)
  sigma      <- sigma(km, h_hat)

  result <- matrix(rep(0, 5 * 2), 5)

  for (i in 1:5) {
    result[i, ] <- statistic(km, sigma, integral, C[[i ]])
  }

  maineffectA       <- result[1,]
  maineffectB       <- result[2,]
  simpleeffectA     <- result[3,]
  simpleeffectB     <- result[4,]
  interactioneffect <- result[5,]


  output <- list(maineffectA = maineffectA, maineffectB = maineffectB, simpleeffectA = simpleeffectA,
               simpleeffectB = simpleeffectB, interactioneffect=interactioneffect, cells = cells)

  names(output[[1]])<-c("statistic", "p_value")
  names(output[[2]])<-c("statistic", "p_value")
  names(output[[3]])<-c("statistic", "p_value")
  names(output[[4]])<-c("statistic", "p_value")
  names(output[[5]])<-c("statistic", "p_value")
  names<-NULL

  names(output[[6]]) <- names

  class(output)<-"nparsurv"

  return (output)
}




