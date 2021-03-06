% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redefine_utils.R
\name{redefine_masks_mask}
\alias{redefine_masks_mask}
\title{IFC_masks Mask Redefinition}
\usage{
redefine_masks_mask(
  masks,
  images,
  to_match_mask = NULL,
  to_replace_mask = NULL,
  ...
)
}
\arguments{
\item{masks}{an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.}

\item{images}{a data.frame containing images definition. Default is missing.}

\item{to_match_mask}{a string with a mask name to use for matching 'masks' names. Default is NULL}

\item{to_replace_mask}{a string of mask name to use for replacing 'masks' names. Default is NULL}

\item{...}{Other arguments to be passed.}
}
\value{
a vector of masks definition or names depending on definition.
}
\description{
Helper to rename a mask within IFC_masks.
}
\keyword{internal}
