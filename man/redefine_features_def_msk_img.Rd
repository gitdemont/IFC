% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redefine_utils.R
\name{redefine_features_def_msk_img}
\alias{redefine_features_def_msk_img}
\title{IFC_features_def Mask or Image Redefinition}
\usage{
redefine_features_def_msk_img(features_def, masks, images, ...)
}
\arguments{
\item{features_def}{an `IFC_features_def` object, or a list containing features definition. Default is missing.}

\item{masks}{an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.}

\item{images}{a data.frame containing images definition. Default is missing.}

\item{...}{Other arguments to be passed.}
}
\value{
a list whose members are:\cr
-features_def, an `IFC_features_def` object, or a list containing features definition\cr
-masks, an `IFC_masks` object or a data.frame containing masks definition and name.\cr
-images, a data.frame containing images definition.
}
\description{
Helper modify features_def according to masks or images redefinition
}
\keyword{internal}
