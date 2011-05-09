##' A Planar Straight Line Graph (PSLG) is a collection of vertices
##' and segments. Segments are edges whose endpoints are vertices in
##' the PSLG, and whose presence in any mesh generated from the PSLG
##' is enforced.
##'
##' @title Create a Planar Straight Line Graph object
##' @param V A 2-column matrix of x-y co-ordinates of vertices. There
##' is one row per vertex.
##' @param VA Matrix of \emph{attributes} which are typically
##' floating-point values of physical quantities (such as mass or
##' conductivity) associated with the nodes of a finite element
##' mesh. When triangulating using \code{\link{triangulate}} these are
##' copied unchanged to the output mesh.
##' @param VB Vector of \emph{boundary markers} of vertices. For each
##' vertex this is 1 if the point should be on a boundary of any mesh
##' generated from the PSLG and 0 otherwise. There should be as many
##' elements in \code{VB} as there are vertices in \code{V}.
##' @param S A 2-column matrix of \emph{segments}. Segments are edges
##' whose endpoints are vertices in the PSLG, and whose presence in
##' any mesh generated from the PSLG is enforced. Each segment refers
##' to the indices in \code{V} of the endpoints of the segment.  There
##' is one row per segment.
##' @param SB Vector of boundary markers of segments. For each segment
##' this is 1 if the segment should be on a boundary of any mesh
##' generated from the PSLG and 0 otherwise. There should be as many
##' elements in \code{SB} as there are segments in \code{S}.
##' @param H 2-column matrix of \emph{holes},  with one hole per
##' row.Holes are specified by identifying a point inside each
##' hole. After the triangulation is formed, Triangle creates holes by
##' eating triangles, spreading out from each hole point until its
##' progress is blocked by PSLG segments; you must be careful to
##' enclose each hole in segments, or your whole triangulation might
##' be eaten away. If the two triangles abutting a segment are eaten,
##' the segment itself is also eaten. Do not place a hole directly on
##' a segment; if you do, Triangle will choose one side of the segment
##' arbitrarily.
##' @return An object containing the input of type \code{pslg} that
##' contains the information supplied in the inputs. This function
##' does some sanity checking of its inputs.
##' @author David Sterratt
pslg <- function(V, VA=NA, VB=NA, S=NA, SB=NA, H=NA) {
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
  
  check.na.nan(V)

  ## Deal with V
  if (ncol(V) != 2) {
    stop("Matrix of verticies V should have 2 columns")
  }

  ## Check that there are no duplicate rows in V
  if (anyDuplicated(V)) {
    stop("Duplicated vertices in V.")
  }

  ## If attributes not specified, set them to zero
  if (any(is.na(VA))) {
    VA <- rep(0, nrow(V))
  }
  
  ## If boudary vertices not specified, set them to 1
  if (is.na(VB)) {
    VB <- rep(1, nrow(V))
  }

  ## Deal with S
  if (any(is.na(S))) {
    S <- cbind(1:nrow(V), c(2:nrow(V),1))
  } else {
    if (ncol(S) != 2) {
    stop("Matrix of segments S should have 2 columns")
    }
  }

  ## If boundary segments not specified, set them to 1
  if (any(is.na(SB))) {
    SB <- rep(1, nrow(S))
  }

  storage.mode(S) <- "integer"
  
  ## Assemble components
  ret <- list(V=V, VA=VA, VB=as.integer(VB),
              S=S, SB=as.integer(SB), H=H)
  class(ret) <- "pslg"
  return(ret)
}

##' Read a Planar Straight Line Graph from a \code{.poly} file, as
##' used in Schwchuk's Triangle library
##' (\url{http://www.cs.cmu.edu/~quake/triangle.poly.html}).
##'
##' @title Read a Planar Straight Line Graph from file
##' @param file File name of \code{.poly} file to read.
##' @return \code{pslg} object. See \code{\link{pslg}}.
##' @author David Sterratt
read.pslg <- function(file) {
  ##   * First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
  ##   * Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
  ##   * One line: <# of segments> <# of boundary markers (0 or 1)>
  ##   * Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
  ##   * One line: <# of holes>
  ##   * Following lines: <hole #> <x> <y>
  ##   * Optional line: <# of regional attributes and/or area constraints>
  ##   * Optional following lines: <region #> <x> <y> <attribute> <maximum area>

  dat <- scan(file, quiet=TRUE)

  ## Read in numbers of vertices, attributes and boundary points
  N.vert <- dat[1]
  N.dims <- dat[2]
  N.attr <- dat[3]
  N.boun <- dat[4]

  ## Read in vertices, attributes and vertex boundaries
  offset <- 4
  line.length <- 3 + N.attr + N.boun

  V <- matrix(NA, N.vert, 2)
  if (N.attr >= 1) {
    VA <- matrix(NA, N.vert, N.attr)
  }
  if (N.boun >= 1) {
    VB <- matrix(NA, N.vert, N.boun)
  } else {
    VB <- NA
  }

  for (i in (1:N.vert)) {
    V[i,] <- dat[offset+((i-1)*line.length)+(2:3)]
    if (N.attr >= 1) {
      VA[i,] <- dat[offset+((i-1)*line.length)+3+(1:N.attr)]
    }
    if (N.boun >= 1) {
      VB[i,] <- dat[offset+((i-1)*line.length)+3+N.attra+(1:N.boun)]
    }
  }

  ## Read in numbers of segments
  offset <- offset + line.length*N.vert
  N.seg <- dat[offset+1]
  N.boun <- dat[offset+2]

  ## Read in segments
  offset <- offset + 2
  line.length <- 3 + N.boun

  S <- matrix(NA, N.seg, 2)
  for (i in (1:N.seg)) {
    S[i,] <- dat[offset+((i-1)*line.length)+(2:3)]
  }

  ## Read in number of holes
  offset <- offset + line.length*N.seg
  N.hole <- dat[offset+1]

  ## Read in holes
  offset <- offset + 1
  H <- matrix(NA, N.hole, 2)
  for (i in (1:N.hole)) {
    H[i,] <- dat[offset+((i-1)*2)+(2:3)]
  }

  return(pslg(V=V, VA=VA, VB=VB, S=S, H=H))
}

##' Plot \code{\link{pslg}} object
##'
##' @title Plot \code{\link{pslg}} object
##' @param p \code{\link{pslg}} object
##' @author David Sterratt
plot.pslg <- function(p) {
  with(p, {
    plot(V)
    segments(V[S[,1],1], V[S[,1],2],
             V[S[,2],1], V[S[,2],2])
  })
}
##' Plots a triangulation object produced with \code{\link{trianglate}}
##'
##' @title Plot a triangulation object produced with \code{\link{triangulate}} 
##' @param t Triangulation object produced with \code{\link{triangulate}}.
##' @author David Sterratt
plot.triangulation <- function(t) {
  with(t, {
    plot(rbind(V))
    segments(V[E[,1],1], V[E[,1],2],
             V[E[,2],1], V[E[,2],2])
    segments(V[S[,1],1], V[S[,1],2],
             V[S[,2],1], V[S[,2],2], col="red")
  })
}

##' Triangulate a planar straight line graph using the Triangle
##' library (\url{http://www.cs.cmu.edu/~quake/triangle.html}).  The
##' triangulation is a constrained conforming Delaunay triangulation
##' in which additional vertices, called Steiner points, can be
##' inserted into segments to improved the quality of the
##' triangulation.  To prevent the insertion of Steiner points on
##' boundary segments, specify \code{Y=1}. If the maximum triangle
##' area \code{a} is specified, the area of each triangle is not
##' allowed to exceed this value. If the the minimum angle \code{q} is
##' specified, no triangle angle is allowed to be below this value.
##'
##' @title Triangulate a Planar Straight Line Graph
##' @param p Planar straight line graph object; see
##' \code{\link{pslg}}.
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
##' @return A object with class \code{triangulation}. This contains
##' the information in the input PSLG, \code{p}, an also contains:
##' \describe{
##' \item{T}{Triangulation specified as 3 column matrix
##' in which each row contains indices in \code{P} of vertices.}
##' \item{E}{Set of edges in the triangulation.}
##' \item{EB}{Boundary markers of edges. For each edge this is 1 if
##' the point is on a boundary of the triangulation and 0
##' otherwise.}}
##' @examples
##' ## Plot PSLG
##' data("A", package="Triangle")
##' plot(A)
##' ## Triangulate the PSLG
##' TA <- triangulate(A)
##' plot(TA)
##' ## Triangulate the PSLG with triangles in which no angle
##' ## is smaller than 20 degrees
##' TA <- triangulate(A, q=20)
##' plot(TA)
##' ## Triangulate the PSLG with triangles in which no triangle has 
##' ## area greater than 0.001
##' TA <- triangulate(A, a=0.001)
##' plot(TA)
##' @author David Sterratt
triangulate <- function(p, a=NULL, q=NULL, Y=FALSE, j=FALSE,
                        V=0, Q=TRUE) {
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
  
  check.na.nan(a)
  check.na.nan(q)
  check.na.nan(V)
  
  ## Call the main routine
  out <- .Call("R_triangulate",
               t(p$V),
               p$VB,
               t(p$S),
               p$SB,
               p$H,
               a,
               q,
               Y,
               j,
               as.integer(V),
               Q)
  names(out) <- c("V", "VB", "T", "S", "SB", "E", "EB")
  class(out) <- "triangulation"
  return(out)
}

# LocalWords:  param Sterratt
