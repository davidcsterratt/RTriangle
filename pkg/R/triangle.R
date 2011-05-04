##' Generate exact Delaunay triangulations, constrained Delaunay
##' triangulations, and high-quality triangular meshes using the
##' Triangle library
##' (\url{http://www.cs.cmu.edu/~quake/triangle.html})
##'
##' Triangulate a planar straight line graph comprising points
##' \code{P} and, optionally, segments \code{S}. The
##' triangulation is a constrained conforming Delaunay triangulation
##' in which additional vertices, called Steiner points, can be
##' inserted into segments to improved the quality of the
##' triangulation.  To prevent the insertion of Steiner points on
##' boundary segments, specify \code{Y=1}. If the maximum triangle area
##' \code{a} is specified, the area of each triangle is not allowed to
##' exceed this value. If the the minimum angle \code{q} is
##' specified, no triangle angle is allowed to be below this value.
##' @title Triangulate a Planar Straight Line Graph
##' @param P Matrix of x-y co-ordinates of vertices, arranged either
##' in rows or columns.
##' @param S Matrix of segments, arranged either in rows or
##' columns. Each segment refers to the indices in \code{P} of the
##' endpoints of the segment.
##' @param a Maximum triangle area. If specified, triangles cannot be
##' larger than this area.
##' @param q Minimum triangle angle in degrees.
##' @param Y If \code{TRUE} prohibits the insertion of Steiner points
##' on the mesh boundary.
##' @param j If \code{TRUE} Jettisons vertices that are not part of
##' the final triangulation from the output.
##' @param V Verbosity level. Specify higher values  for more detailed
##' information about what the Triangle library is doing.
##' @param Q If \code{TRUE} suppresses all explanation of what the
##' Triangle library is doing, unless an error occurs. 
##' @return \item{P}{Set of vertices in the triangulation.}
##' \item{PB}{Boundary markers of vertices. For each vertex this is 1
##' if the point is on a boundary of the triangulation and zero
##' otherwise.}
##' \item{T}{Triangulation specified as 3 column matrix
##' in which each row contains indices in \code{P} of vertices.}
##' \item{S}{Set of segments (enforced edges) in the triangulation.}
##' \item{SB}{Boundary markers of segments. For each segment this is 1
##' if the point is on a boundary of the triangulation and 0
##' otherwise.}
##' \item{E}{Set of edges in the triangulation.}
##' \item{EB}{Boundary markers of edges. For each edge this is 1 if
##' the point is on a boundary of the triangulation and 0
##' otherwise.}
##' @author David Sterratt
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

# LocalWords:  param Sterratt
