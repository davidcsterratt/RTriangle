triangulate <- function(P, S=NULL, a=NULL, q=NULL, Y=FALSE, j=FALSE,
                        V=0, Q=FALSE) {
  ## It is necessary to check for NAs and NaNs, as the triangulate C
  ## code crashes if fed with them
  check.na.nan <- function(x) {
    if (!is.null(x)) {
      if (any(is.nan(x))) {
        stop(paste("NaN in", deparse(substitute(x))))
      }
      if (any(is.na(x))) {
        stop(paste("NA in", deparse(substitute(x))))
      }
    }
  }
  
  check.na.nan(P)
  check.na.nan(S)
  check.na.nan(a)
  check.na.nan(q)
  check.na.nan(V)

  ## Deal with P
  if (ncol(P) == 2) {
    P <- t(P)
  }
  ## Check that there are no duplicate rows in P
  if (anyDuplicated(t(P))) {
    stop("Duplicated points in P.")
  }
  PB <- rep(1, ncol(P))

  ## Deal with S
  if (is.null(S)) {
    S <- rbind(1:ncol(P), c(2:ncol(P),1))
  } else {
    if (ncol(S) == 2) {
      S <- t(S)
    }
  }
  SB <- rep(1, ncol(S))
  
  ## Call the main routine
  out <- .Call("R_triangulate",
               P,
               as.integer(PB),
               as.integer(S),
               as.integer(SB),
               a,
               q,
               Y,
               j,
               as.integer(V),
               Q)
  names(out) <- c("P", "PB", "T", "S", "SB", "E", "EB")
  return(out)
}
