\name{pslg}
\alias{pslg}
\title{Create a Planar Straight Line Graph object}
\usage{
pslg(P, PB = NA, PA = NA, S = NA, SB = NA, H = NA)
}
\arguments{
  \item{P}{A 2-column matrix of x-y co-ordinates of
  vertices. There is one row per vertex.}

  \item{PB}{Vector of \emph{boundary markers} of vertices.
  For each vertex this is 1 if the point should be on a
  boundary of any mesh generated from the PSLG and 0
  otherwise. There should be as many elements in \code{VB}
  as there are vertices in \code{V}.}

  \item{PA}{Matrix of \emph{attributes} which are typically
  floating-point values of physical quantities (such as
  mass or conductivity) associated with the nodes of a
  finite element mesh. When triangulating using
  \code{\link{triangulate}} these are copied unchanged to
  existing points in the output mesh and each new Steiner
  point added to the mesh will have quantities assigned to
  it by linear interpolation.}

  \item{S}{A 2-column matrix of \emph{segments} in which
  each row is a \emph{segment}. Segments are edges whose
  endpoints are vertices in the PSLG, and whose presence in
  any mesh generated from the PSLG is enforced. Each
  segment refers to the indices in \code{V} of the
  endpoints of the segment. By default the segments are not
  specified (\code{NA}), in which case the convex hull of
  the vertices are taken to be the segments. Any vertices
  outside the region enclosed by the segments are eaten
  away by the triangulation algorithm. If the segments do
  not enclose a region the whole triangulation may be eaten
  away.}

  \item{SB}{Vector of boundary markers of segments. For
  each segment this is 1 if the segment should be on a
  boundary of any mesh generated from the PSLG and 0
  otherwise. There should be as many elements in \code{SB}
  as there are segments in \code{S}.}

  \item{H}{2-column matrix of \emph{holes}, with one hole
  per row.Holes are specified by identifying a point inside
  each hole. After the triangulation is formed, Triangle
  creates holes by eating triangles, spreading out from
  each hole point until its progress is blocked by PSLG
  segments; you must be careful to enclose each hole in
  segments, or your whole triangulation might be eaten
  away. If the two triangles abutting a segment are eaten,
  the segment itself is also eaten. Do not place a hole
  directly on a segment; if you do, Triangle will choose
  one side of the segment arbitrarily.}
}
\value{
An object containing the input of type \code{pslg} that
contains the information supplied in the inputs. This
function does some sanity checking of its inputs.
}
\description{
A Planar Straight Line Graph (PSLG) is a collection of
vertices and segments. Segments are edges whose endpoints
are vertices in the PSLG, and whose presence in any mesh
generated from the PSLG is enforced.
}
\author{
David Sterratt
}

