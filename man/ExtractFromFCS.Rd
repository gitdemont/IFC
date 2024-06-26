% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs.R
\name{ExtractFromFCS}
\alias{ExtractFromFCS}
\title{FCS File Reader}
\source{
Data File Standard for Flow Cytometry, version FCS 3.2 from Spidlen J. et al. available at \doi{10.1002/cyto.a.24225}.
}
\usage{
ExtractFromFCS(fileName, ...)
}
\arguments{
\item{fileName}{path(s) of file(s). If multiple files are provided they will be merged and 
populations will be created to identify each single file within returned `IFC_data` object.}

\item{...}{other arguments to be passed to readFCS function, with the exception of 'options$text_only'.\cr
Experimental (as of v0.2.1.300), ExtractFromFCS could handle FCS 3.2 files. However, it is important to note that R has no native support for 64bits unsigned integers which are defined in the FCS 3.2 standard.
So, large integers are extracted as double (8 bytes) and precision loss will happen for > 2^53 integers on 64bits platforms.}
}
\value{
A named list of class `IFC_data`, whose members are:\cr
-description, a list of descriptive information,\cr
-Merged_fcs, character vector of path of files used to create fcs, if input was a merged,\cr
-Keywords, a named-list of keywords values, only keywords from 1st 'fcs' segment will be retrieved\cr
-fileName, path of fileName input,\cr
-fileName_image, path of .cif image fileName is referring to,\cr
-features, a data.frame of features,\cr
-features_def, a describing how features are defined,\cr
-graphs, a list of graphical elements found,\cr
-pops, a list describing populations found,\cr
-regions, a list describing how regions are defined,\cr
-images, a data.frame describing information about images,\cr
-offsets, an integer vector of images and masks IFDs offsets,\cr
-stats, a data.frame describing populations count and percentage to parent and total population,\cr
-checksum, a checksum integer.
}
\description{
Extracts data from Flow Cytometry Standard (FCS) Files.
}
