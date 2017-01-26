################################
# Contrast Matrices            #
################################

#==============================================================================================================

# input:  - a,b = number of levels of Variable A respectively B
#         - ind = logical vector of length 5 which matrix should be computed

# output: list with the wished contrast matrices as in AKritas and Brunner

#==============================================================================================================

contrast <- function(a, b, ind) {

  contrast <- list(NULL)
  c    <- a*b
  vec  <- c(1:c)
  A    <- matrix(vec, a, b, byrow=TRUE)
  perm <- c(A)

  # Contrast matrix for main effect of A
  if (ind[1] == 1) {
  M_a      <- cbind(rep(1, a - 1), - diag(a - 1))
  C_A      <- kronecker(M_a, (1 / b) * t(rep(1, b)))
  C_A      <- list(C_A = C_A)
  contrast <- c(contrast, C_A)
  }

  # Contrast matrix for main effect of B
  if (ind[2] == 1) {
  C_B      <- t(matrix(contrast$C_A[, perm], c, byrow=TRUE))
  contrast <- c(contrast, list(C_B= C_B))
  }

  # Contrast matrix for no simple effects of A
  if (ind[3] == 1) {
  M_a      <- cbind(rep(1, a - 1), - diag(a - 1))
  C_A_B    <- kronecker(M_a, diag(b))
  C_A_B    <- list(C_A_B = C_A_B)
  contrast <- c(contrast, C_A_B)
  }

  # Contrast matrix for no simple effects of B
  if (ind[4] == 1) {
  C_B_A    <- t(matrix(contrast$C_A_B[, perm], c, byrow=TRUE))
  contrast <- c(contrast, C_B_A = list(C_B_A))
  }

  # Contrast matrix for interaction effects
  if (ind[5] == 1) {
    M_a      <- cbind(rep(1, a - 1), - diag(a - 1))
    M_b      <- cbind(rep(1, b - 1), - diag(b - 1))
    C_AB     <- kronecker(M_a, M_b)
    C_AB     <- list(C_AB = C_AB)
    contrast <- c(contrast, C_AB)
  }

  return(contrast[-1])
}
