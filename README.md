# Tools for Imaging Flow Cytometry (IFC)

## INSTALLATION

### On Windows (tested 7 and 10)

- install [Rtools for windows](https://cran.r-project.org/bin/windows/Rtools/)

- ensure that Rtools compiler  in present in Windows PATH

In R console, if you installed Rtools directly on C: you should see something like C:\\Rtools\\bin and C:\\Rtools\\mingw_32\\bin by typing: unlist(strsplit(Sys.getenv("PATH"), ";"))

Otherwise, try shell('setx PATH "C:\\Rtools\\bin"') and shell('setx PATH "C:\\Rtools\\mingw_32\\bin"')

You may need to reboot to operate changes

### On MacOS (tested on High Sierra)

- IFC package seems to installed without additional requirement.

### On Linux (tested on Ubuntu 16.04 LTS)

- IFC package is dependent on tiff package. It requires tiff libraries which may not be on the system.

sudo apt-get install libtiff-dev 

### in R

- install R dependencies require for IFC package : "Rcpp", "RcppProgress", "xml2", "png", "tiff", "jpeg", "utils", "grid", "gridExtra", "lattice", "latticeExtra", "KernSmooth", "DT", "visNetwork"

install.packages(c("Rcpp", "RcppProgress", "xml2", "png", "tiff", "jpeg", "utils", "grid", "gridExtra", "lattice", "latticeExtra", "KernSmooth", "DT", "visNetwork"))

- install "devtools", to install IFC package from github devtools is needed.

install.packages("devtools")

- install IFC

library(devtools)

install_github(repo = "gitdemont/IFC", ref = "master", dependencies = FALSE)


## USAGE

library(IFC)

help(IFC)

## DETAILS

- parse files .daf

use: **readIFC()** or **ExtractFromDAF()**

This allows retrieving several information from files like: 

> features defined and their values, graphs from the analysis worksheet,

> masks defined,

> regions drawn,

> populations created.

As a consequence ones can get for each cell the population it belongs to. 
One main advantage is that you don't have to use export / extract .cif or .fcs or .txt manually for each subpopulation you are interested in.
Having features values for each cells and knowing which populations a cell is belonging to, one can pass it to a supervised machine learning (ML) algorithm (they are many in R) to create a model fitting the data (IDEAS only provides Rd Ratio).

- parse files .rif, .cif

use: **readIFC()** or **ExtractFromXIF()**

This allows to directly extract raw / compensated images and associated masks for all acquired channels. These images can be then easily exported to R arrays or TIFF / JPEG / PNG images.
For instance, we can programmatically export images or images values from a desired population indentified in a daf file.
Accessing images allow image export to create beautiful gallery export but also to pass it to deep learning frameworks like Tensorflow or MXNet (e.g. Convolutional Neural Network)

- write in .daf

use: **writeIFC()** or **ExportToDAF()**

Once data are parsed and treated new elements can be injected in daf file, e.g. a ML model is trained on a daf file and used to predict populations in other daf files.
Then these predicted populations can be injected so as to be checked using IDEAS in addition to their predicted probability values (or other features like PCA/t-SNE dimensions).

- subset or merge .rif, .cif,

use: **writeIFC()** or **ExportToXIF()**
Once data are parsed and treated new elements can be used to subset or merge raw / compensated images files

- Scale up productivity; here are several (among others) functions that have been created to:

> display cells / gallery of cells in R, use **DisplayGallery()**,

> export gallery of cells (with desired channels) to image files (tiff, png, jpeg, …), use **ExportToGallery()**,

> export gallery of cells (with desired channels) to Numpy arrays, use **ExportToNumpy()**,

> export graphs and associated stats from analysis worksheet of files, use **ExportToReport()**,

> view in R / export population hierarchy, use **popsNetwork()**,

> create batch demand to be processed by IDEAS, use **ExportToBATCH()**.

- Other low level functions are very useful see ?IFC, for how to use them.

**getInfo()** for retrieving rich information about file,

**objectExtract()** for extracting individual image or mask within .rif or .cif file,

**plotGraph()** and  **autoplot()** for displaying graphical elements of analysis worksheet.

## DISCLAMER

- You are using this package **on your own risk!**

- We do not guarantee **privacy nor confidentiality**.

- This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**. In no event shall the copyright holders or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.
