% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_add_pop_sample.R
\name{data_add_pop_sample}
\alias{data_add_pop_sample}
\title{Population Sampling}
\usage{
data_add_pop_sample(obj, pop, size, new_name, random_seed = NULL, ...)
}
\arguments{
\item{obj}{an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).}

\item{pop}{name of the population to sample.}

\item{size}{a non-negative integer giving the number of items to choose.}

\item{new_name}{name of the exported population.}

\item{random_seed}{a list of elements to pass to \link[base]{set.seed} or a single value, interpreted as an integer, or NULL to be used when 'add_noise' is set to TRUE. Default is NULL.
Note that NA_integer_ or list(seed = NA_integer_) can be used to not call \link[base]{set.seed} at all.}

\item{...}{Other arguments to be passed.}
}
\value{
an IFC_data object with sampled pop added.
}
\description{
Creates a sample from a population
}
\details{
population is exported as tagged population.
}
\keyword{internal}
