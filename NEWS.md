# NEWS
## 0.1.2
- allow spatial correction directly from objectExtract thanks to the addition of spatial_correction parameter in objectParam

- new function to convert XIF to TIFF (+ underlying hpp/cpp functions)

- remove TILE support and only use STRIP in mergeXIF and susbsetXIF since it seems to never been used

- remove hard dependency on object number retrieved from first IFD in getIFD
allowing to get IFDs from other files than RIF / CIF (i.e. TIFF)

- make decomp.hpp more robust to invalid nbytes, imgWidth, imgHeight, nb_channels input 
fix potential off-by-one issue that could cause buffer overrun in code although it should never be reached 

- new function to deal with compensation

- new function to apply spatial offsets correction on images

- new function to retrieve ASSIST tests values

- fix bug in setting Raw Number features definition when using ExtractFromXIF

- add new color in allowed colors palette "Control" now matches "Gray81"

- fix bug that disabled styling edition of population "All"

- add maxpoints parameter in buildGraph to control the number of cells displayed in 2D graphs

- create read/write/apply GatingStrategy family functions to save/retrieve gating strategy
(association of regions / populations and graphs)

- better compute graph dependency in data_rm_ family functions

- add adjust_graph parameter to control if graph should be removed or if it should be modified when a

- improve getAborted to retrieve file path of aborted elements in batch

- fix bug in shiny progress bar which did not initialize with title / detail

- improve mergeXIF and subsetXIF functions (speed gain + no more dependency on seek for reading current position)

## 0.1.1
- CRAN release

## 0.1.0
- fix bug when using verbose in Extract Images/Masks + ExportTo/Display Gallery/Numpy family functions

- fix bug in objectDisplay when image has more than one class

- fix bug in ExportToNumpy object size was not correctly set

- pass all cpp functions to header to allow more portability (side effect package installation is faster)

- reintegrate corrected version of num_to_string function using format() did not produce desired results

- add functions family to remove features, regions and / or populations from IFC_data object

- allow to create graph report even if one graph produces an error
also allow to produce graph with base population with no event (change stop to warning)
also allow to produce density graph with aberrant shown pop filled (change stop to warning)

## 0.0.9
- allow graph conversion from lattice to base

- modify plotGraph return object to allow to pass it to base plot + better handle graphics device

- fix solaris installation error