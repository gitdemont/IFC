# NEWS
## 0.1.3
- don't send error on partial mask retrieval, flag removal with "invalid" instead

## 0.1.2
- CRAN release

- add color-gradient based on 3rd feature in density plot + add this 3rd feature as subtitle

- fix issue with axis tick labels on Solaris (polygon edge not found (zero-width or zero-height?))

- fix issue in internal addText() and hpp_mark() when added text does not fit into image.

- allow stats export with data_to_DAF and daf creation from rif file reading

- make changes to num_to_string to optimize numeric to string conversion when exporting file

- fix bug with rif reading (ExtractFromXIF) that was also extracting not acquired channels

- fix images values retrieval in getImagesValues (only one channel was retrieved instead of all)

- handle "0" color value

- retrieve "EVmode" information with getInfo

- add experimental FCS reader / writer

- improve axes labels positioning

- modify how transformation instruction is parsed to allow for more transformations in the future

- apply STRICT_R_HEADERS patch for Rcpp and change PI to M_PI

- fix bug with base64id that was not quoted which could lead id to be truncated on space characters

- modify string_utils to better escape and recognize regions, populations and features names, particularly when they contain | or when names are constituted of repeated patterns

- add redefine_ family internal functions to allow images, masks, and features renaming and their propagation across upstream dependency

- add support to more XIF files thanks to new internal testXIF function that allow to identify files with no mask

- allow spatial correction directly from objectExtract thanks to the addition of spatial_correction parameter in objectParam

- new experimental internal function to convert XIF to uncompressed TIFF (+ underlying hpp/cpp functions)

- remove TILE support and only use STRIP in mergeXIF and susbsetXIF since it seems to never been used

- remove hard dependency on object number retrieved from first IFD in getIFD allowing to get IFDs from other files than RIF / CIF (i.e. TIFF)

- make decomp.hpp more robust to invalid nbytes, imgWidth, imgHeight, nb_channels input fix potential off-by-one issue that could cause buffer overrun in code although it should never be reached 

- new function to apply spatial offsets correction on images

- new function to retrieve ASSIST tests values

- fix bug in setting Raw Number features definition when using ExtractFromXIF

- add new color in allowed colors palette "Control" now matches "Gray81"

- fix bug that disabled styling edition of population "All"

- add maxpoints parameter in buildGraph to control the number of cells displayed in 2D graphs

- create experimental read/write/apply GatingStrategy family functions to save/retrieve gating strategy (association of regions / populations and graphs) with partial support for GatingML

- better compute graph dependency in data_rm_ family functions

- add adjust_graph parameter to control if graph should be removed or if it should be modified if possible

- improve getAborted to retrieve file path of aborted elements in batch

- fix bug in shiny progress bar which did not initialize with title / detail

- improve mergeXIF and subsetXIF functions (speed gain + no more dependency on seek for reading current position)

#### This leads to the following visible changes for the user
*data_to_DAF now allows to create daf file with statistics table*

*data_to_DAF now allows to create daf file from an `IFC_data` object extracted from a rif file*

*ExtractFromXIF does not extract channels that were not acquired anymore*

*getInfo returned object gains an additional `evmode` value*

*buildGraph gains a new `maxpoints` parameter*

*data_rm_ functions gain a new `adjust_graph` parameter*

*subsetOffsets and getOffsets returned object gain an additional `test` attribute*

*getInfo now returns `XIF_test` value in addition to former elements*

*objectParam gains a new `spatial_correction` parameter which when TRUE adds a 2 new columns, namely `spatial_X` and `spatial_Y`, to the returned data.frame in 'channels'*

*paletteIFC now maps Control to Gray81*

*axes display have been slightly changed to show labels at `pretty` positions when LinLog transformation is used*

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
