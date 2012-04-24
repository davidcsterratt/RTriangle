\name{triangulate}
\alias{triangulate}
\title{Triangulate a Planar Straight Line Graph}
\usage{
  triangulate(p, a = NULL, q = NULL, Y = FALSE, j = FALSE,
    D = FALSE, S = Inf, V = 0, Q = TRUE)
}
\arguments{
  \item{p}{Planar straight line graph object; see
  \code{\link{pslg}}.}

  \item{a}{Maximum triangle area. If specified, triangles
  cannot be larger than this area.}

  \item{q}{Minimum triangle angle in degrees.}

  \item{Y}{If \code{TRUE} prohibits the insertion of
  Steiner points on the mesh boundary.}

  \item{j}{If \code{TRUE} jettisons vertices that are not
  part of the final triangulation from the output.}

  \item{D}{If \code{TRUE} produce a conforming Delaunay
  triangulation. This ensures that all the triangles in the
  mesh are truly Delaunay, and not merely constrained
  Delaunay.  This option invokes Ruppert's original
  algorithm, which splits every subsegment whose diametral
  circle is encroached.  It usually increases the number of
  vertices and triangles.}

  \item{S}{Specifies the maximum number of added Steiner
  points.}

  \item{V}{Verbosity level. Specify higher values for more
  detailed information about what the Triangle library is
  doing.}

  \item{Q}{If \code{TRUE} suppresses all explanation of
  what the Triangle library is doing, unless an error
  occurs.}
}
\value{
  A object with class \code{triangulation}. This contains
  the information in the same format as the PSLG, \code{p},
  with an updated list of points \code{P} and point
  attributes \code{PA}, along with the following variables:

  \item{\code{T}}{Triangulation specified as 3 column
  matrix in which each row contains indices in \code{P} of
  vertices.} \item{\code{E}}{Set of edges in the
  triangulation.} \item{\code{EB}}{Boundary markers of
  edges. For each edge this is 1 if the point is on a
  boundary of the triangulation and 0 otherwise.}
  \item{\code{VP}}{The points of the Voronoi tessalation as
  a 2-column matrix} \item{\code{VE}}{Set of edges of the
  Voronoi tessalation. An index of -1 indicates an infinite
  ray.} \item{\code{VN}}{Directions of infinite rays of
  Voroni tessalation as a 2-column matrix with the same
  number of rows as \code{VP}.} \item{\code{VA}}{Matrix of
  \emph{attributes} associated with the polygons of the
  Voronoi tessalation.}
}
\description{
  Triangulate a planar straight line graph using the
  Triangle library
  (\url{http://www.cs.cmu.edu/~quake/triangle.html}).  The
  triangulation is a constrained conforming Delaunay
  triangulation in which additional vertices, called
  Steiner points, can be inserted into segments to improved
  the quality of the triangulation.  To prevent the
  insertion of Steiner points on boundary segments, specify
  \code{Y=1}. If the maximum triangle area \code{a} is
  specified, the area of each triangle is not allowed to
  exceed this value. If the the minimum angle \code{q} is
  specified, no triangle angle is allowed to be below this
  value.
}
\examples{
## Create an object with a concavity
p <- pslg(P=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
          S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)))
## Plot it
plot(p)
## Triangulate it
tp <- triangulate(p)
plot(tp)
## Triangulate it subject to minimum area constraint
tp <- triangulate(p, a=0.01)
plot(tp)
## Load a data set containing a hole
A <- read.pslg(file.path(system.file(package = "RTriangle"), "extdata", "A.poly"))
plot(A)
## Produce a constrained Delaunay triangulation of the PSLG
tA <- triangulate(A, Y=TRUE)
plot(tA)
## Produce a conforming Delaunay triangulation of the PSLG
tA <- triangulate(A, D=TRUE)
plot(tA)
## Triangulate the PSLG with triangles in which no angle
## is smaller than 20 degrees
tA <- triangulate(A, q=20)
plot(tA)
## Triangulate the PSLG with triangles in which no triangle has
## area greater than 0.001
tA <- triangulate(A, a=0.001)
plot(tA)
}
\author{
  David Sterratt
}

