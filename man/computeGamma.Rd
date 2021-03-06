% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma.R
\name{computeGamma}
\alias{computeGamma}
\title{Image Gamma Computation}
\usage{
computeGamma(V)
}
\arguments{
\item{V}{channel display properties containing 'xmin', 'xmax', 'xmid' and 'ymid'.}
}
\description{
Computes image gamma transformation value.
}
\details{
Thanks to 
http://www.libpng.org/pub/png/book/chapter10.html
image_sample = light_out ^ gamma
it is said that
Once again, bear in mind that light_out and image_sample are scaled to the interval between 0 and 1;
that is, if the sample depth is 8 bits, the file samples range between 0 and 255, so image_sample is
obtained by dividing a given file sample by 255, in floating-point arithmetic. 
So,
image_sample = ymid and its range is [0,255]
light_out = xmid and its range is [xmin,xmax]
we have ymid / 255 = ((xmid - xmin) / (xmax - xmin)) ^ gamma
log(ymid / 255) = log((xmid - xmin) / (xmax - xmin)) * gamma
gamma = log(ymid / 255) / log((xmid - xmin) / (xmax - xmin))
}
\keyword{internal}
