<!-- badges: start -->
<!--[![R build status](https://github.com/gitdemont/IFC/workflows/R-CMD-check/badge.svg)](https://github.com/gitdemont/IFC/actions)-->
[![](https://img.shields.io/github/last-commit/gitdemont/IFC.svg)](https://github.com/gitdemont/IFC/commits/master)
[![runiverse IFC status badge](https://gitdemont.r-universe.dev/IFC/badges/version)](https://gitdemont.r-universe.dev/IFC)
[![cran-badge](https://www.r-pkg.org/badges/version/IFC)](https://CRAN.R-project.org/package=IFC)
[![Dependencies](https://tinyverse.netlify.app/badge/IFC)](https://CRAN.R-project.org/package=IFC)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

[![](https://img.shields.io/github/languages/code-size/gitdemont/IFC.svg)](https://github.com/gitdemont/IFC)
[![](http://cranlogs.r-pkg.org/badges/grand-total/IFC?color=blue)](https://cran.r-project.org/package=IFC)
[![](http://cranlogs.r-pkg.org/badges/last-month/IFC?color=green)](https://cran.r-project.org/package=IFC)
[![](http://cranlogs.r-pkg.org/badges/last-week/IFC?color=yellow)](https://cran.r-project.org/package=IFC)
<!-- badges: end -->

# Tools for Imaging Flow Cytometry (IFC)

## INSTALLATION (from CRAN)

The easiest way to install `IFC` package is:

```R
install.packages("IFC")
```

## INSTALLATION (from **in-dev** github master branch)

Current in development branch can be installed from github and requires compilation

#### On Windows (tested 7 and 10)

- install [Rtools for windows](https://cran.r-project.org/bin/windows/Rtools/)

- ensure that Rtools compiler is present in Windows PATH

In R console, if you installed Rtools directly on C: you should see something like C:\\Rtools\\bin and C:\\Rtools\\mingw_32\\bin 

```R
print(unlist(strsplit(Sys.getenv("PATH"), ";")))
```

Otherwise, try to set it

```R
# e.g.
shell('setx PATH "C:\\Rtools\\bin"')
shell('setx PATH "C:\\Rtools\\mingw_32\\bin"')
```

You may need to reboot to operate changes

#### On MacOS (tested on High Sierra)

- IFC package seems to install without additional requirement.

#### On Linux (tested on Ubuntu 16.04 LTS)

- tiff package requires tiff libraries which may not be on the system.

```terminal
sudo apt-get install libtiff-dev
```

#### in R

- install R dependencies required for IFC package
"Rcpp", "xml2", "utils", "grid", "gridExtra", "gridGraphics", "lattice", "latticeExtra", "KernSmooth", "DT", "visNetwork"

```R
install.packages(c("Rcpp", "xml2", "utils", "grid", "gridExtra", "gridGraphics", "lattice", "latticeExtra", "KernSmooth", "DT", "visNetwork"))
```

- install "remotes", to install IFC package from github remotes is needed.

```R
install.packages("remotes")
```

- install IFC

```R
remotes::install_github(repo = "gitdemont/IFC", ref = "master", dependencies = FALSE)
```

## USAGE

Several examples in IFC package are dependent on data files that can be found in dedicated IFCdata package.

To install IFCdata package and run examples in IFC:

```R
install.packages("IFCdata", repos = "https://gitdemont.github.io/IFCdata/", type = "source")
```

## DETAILS

- parse files .daf, .rif, .cif

use: **readIFC()** or **ExtractFromDAF()** or **ExtractFromXIF()**

This allows retrieving several information from files like at `single cell` level: 

`features` defined and their values

`masks` defined and their definitions,

`regions` drawn and their vertices,

`populations` created and know whether a cell belongs to this population or not,

`graphs` from the analysis worksheet.

As a consequence ones can get for each cell the population it belongs to. 
One main advantage is that you don't have to use export / extract .cif or .fcs or .txt manually for each subpopulation you are interested in.
Having features values for each cells and knowing which populations a cell is belonging to, one can pass it to a supervised machine learning (ML) algorithm (they are many in R) to create a model fitting the data.

- write in .daf

use: **writeIFC()** or **ExportToDAF()**

Once data are parsed and treated new elements can be injected in daf file, e.g. a ML model is trained on a daf file and used to predict populations in other daf files.
Then these predicted populations can be injected so as to be checked using IDEAS in addition to their predicted probability values (or other features like PCA/t-SNE/UMAP dimensions).

- subset or merge .rif, .cif,

use: **writeIFC()** or **ExportToXIF()**

Once data are parsed and treated new elements can be used to subset or merge raw / compensated images files

- Scale up productivity; here are several (among others) functions that have been created:

**ExtractFromFCS** or **ExportToFCS**, to read FCS and create FCS files,

**DisplayGallery()**, displays cells / gallery of cells in R,

**ExtractImages_toFile()**, extracts cells (with desired channels) to image files (tiff, png, jpeg, …),

This allows direct extraction of raw / compensated images and associated masks for all acquired channels. These images can be then easily exported to R arrays or TIFF / JPEG / PNG images.
For instance, we can programmatically export images from a desired population identified in a daf file.
Accessing images allows to pass them to deep learning frameworks like tensorflow (e.g. Convolutional Neural Network)

**ExportToReport()**, exports graphs and associated stats from analysis worksheet of files,

**popsNetwork()**, allows visualization in R / exports population hierarchy,

**ExportToBATCH()**, creates batch demand to be processed by IDEAS,

- Other low level functions are very useful see ?IFC, for how to use them.

**getInfo()** for retrieving rich information about file,

**objectExtract()** for extracting individual image or mask within .rif or .cif file,

**plotGraph()** and  **autoplot()** for displaying graphical elements of analysis worksheet,

**data_add_*** and **data_rm_*** functions for data object manipulation,

**data_to_DAF** to export object manipulated in R to .daf file.

## DISCLAIMER

- You are using this package **on your own risk!**

- We do not guarantee **privacy nor confidentiality**.

- This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**. In no event shall the copyright holders or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
