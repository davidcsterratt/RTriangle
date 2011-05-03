##if (is.loaded("R_triangulate")) dyn.unload("R_triangle.so")
##dyn.load("R_triangle.so")

triangulate <- function(P, S=NULL, a=NULL, q=NULL, Y=FALSE, j=FALSE,
                        V=0, Q=FALSE) {
  if (ncol(P) == 2) {
    P <- t(P)
  }
  PB <- rep(1, ncol(P))
  if (is.null(S)) {
    S <- rbind(1:ncol(P), c(2:ncol(P),1))
  } else {
    if (ncol(S) == 2) {
      S <- t(S)
    }
  }
  SB <- rep(1, ncol(S))

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
