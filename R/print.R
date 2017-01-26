###########################
# Print                   #
###########################

#=========================================================
# input: - statistic
#        - p-value
#        - effects
#=========================================================


#'@export

print.nparsurv<- function(x,...) {
  cat("___Tests___\n")
  cat( "\n", "test", "\t", "\t", "\t", "statistic", "\t", "p-value", "\n", "\n")
  cat("Main effect of A","\t", x[[1]][1],"\t", x[[1]][2], "\n")
  cat("Main effect of B", "\t", x[[2]][1],"\t", x[[2]][2], "\n")
  cat("Simple effect of A","\t", x[[3]][1],"\t", x[[3]][2], "\n")
  cat("Simple effect of B","\t", x[[4]][1],"\t", x[[4]][2], "\n")
  cat("Interaction effect","\t", x[[5]][1],"\t", x[[5]][2], "\n", "\n")
}

