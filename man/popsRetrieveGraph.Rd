% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popsRetrieveGraph.R
\name{popsRetrieveGraph}
\alias{popsRetrieveGraph}
\title{Graph Retrieval from Graphical IFC_pops}
\usage{
popsRetrieveGraph(obj, pops, vis2D = "density", all_siblings = FALSE)
}
\arguments{
\item{obj}{an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).}

\item{pops}{names of graphical populations present in 'obj'. Note that they should be siblings.}

\item{vis2D}{when original graph is not an histogram, whether to display it as "scatter" or "density". Default is "density".}

\item{all_siblings}{whether to add all 'pop' siblings in the graph. Default is FALSE.}
}
\value{
a list of parameters needed to build an IFC graph.
}
\description{
Retrieves the graph a graphical population originate from
}
\keyword{internal}
