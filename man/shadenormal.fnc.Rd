\name{shadenormal.fnc}
\alias{shadenormal.fnc}
\title{Shade rejection region for normal probability density function}
\description{
  This function plots the standord normal probability density function
  and shades the rejection region.  
}
\usage{
shadenormal.fnc(qnts = c(0.025, 0.975))
}
\arguments{
  \item{qnts}{A numeric vector with the Z-scores of the boundaries
      of the lower and upper rejection regions.}
}
\value{
  A plot on the graphics device.

	Type \code{shadenormal.fnc} to see the code.  The polygon() function used for
	the shaded areas takes a sequence of X and Y coordinates, connects the
	corresponding points, and fills the area(s) enclosed with a specified color.
	To understand the use of {\tt polygon()}, one can best think of making a
	polygon with a set of pins, a thread, and a board.  Outline the polygon by
	placing the pins on the board at the corners of the polygon. First fasten the
	thread to one of the pins, then connect the thread to the second pin, from
	there to the third pin, and so on, until the first pin has been reached.
	What {\tt polygon()} requires as input is a vector of the X-coordinates of
	the pins, and a vector of their Y-coordinates.  These coordinates should be
	in exactly the order in which the thread is to be connected from pin to pin.  

  For shading the left rejection area, we specify the vectors of X and Y
  coordinates, beginning at the leftmost point of the tail, proceding to the
  right edge of the shaded area, then up, and finally to the left and down to
  the starting point, thereby closing the polygon.  The X-coordinates are
  therefore specified from left to right, and then from right to left.  The
  corresponding Y-coordinates are all the zeros necessary to get from $-3$ to
  $1.96$ (the default, qnorm(0.025)), and then the Y-coordinates of the 
  density in reverse order to return to where we began.  


}
\author{R. H. Baayen}
\examples{
shadenormal.fnc()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
