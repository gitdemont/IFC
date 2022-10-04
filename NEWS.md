# NEWS
## 0.1.6.xxx

- Redefine (Internally)

improve and document logic behind **switch_channel** (used to allow changing all images/masks/features definition of one channel to another one)

create **swap_channel** function (which use **switch_channel**) to operate channel swap (basically, switch to one channel, switch to other channel, and merge)

- Definition Parsing (Internally)

Internally, the framework used for parsing definitions (features/pops/masks) has been consolidated to enhance accuracy and security. This is mostly based on **splitn** (for pops and features) which split definition into chunks according to allowed names and operators and ensure 1-to-1 replacement of them with unique random names created thanks to **gen_altnames** and **random_name**

create dedicated **gen_altnames** function to map allowed names to unique id

adapt **splitn** to handle features definition, simplify logic, make use of **gen_altnames** ancillary function, and avoid unnecessary unique id creation thanks to new `alt_names` argument and add `scalar` and `dsplit` arguments

add dedicated **split_feat** (based on **splitn** and **get_altnames**) for features to avoid repeated splitting loop

- Features (obj$features, in `IFC_data` object)

create **featureIFC** for static features types/names (used to validate definition and also determine name of a feature)

create **gseq**/**cpp_seqmatch** to match string sequence, used by **feature_namer**

modify **feature_namer** (used to determine name of a feature)

"combined" features resulted in wrong values when extracted (dedicated internal **get_feat_value** and **getFeaturesValues** are now here to solve this)

**data_rm_features** caused the removal of more features than it should (due to wrong pattern used for features names matching)

- Populations (obj$pops, in `IFC_data` object)

fix potential bug with NAs / NaNs values, e.g. when graphically defined on features containing NAs / NaNs, the resulted event should be 'FALSE' and not 'NA' to flag that it does **NOT** belong to the population

limit the scope of evaluation for features and populations computation

adding/exporting "tagged" population containing NAs / NaNs value is not allowed

warning/error messages are now more informative during population add/compute/export

- Stats

improve and dry code for stats extraction / computation

group internal **buildStats** and **statsCompute** functions into **buildStats**

add internal functions to extract and generate Statistics Report

fix bug leading to error (pops reported count was NAs) when features contained NAs

fix bug in **plotGraph**, **ExportToReport**, and **BatchReport** in feature summary computation for graphs containing overlay of populations (1D and 2D) for values of the region(s) drawn

- Graphs

catch errors in density computation for histogram and send warning instead but allow graph to be created as far as possible (already catched in **plot_lattice** but not with **plot_base** nor **plot_raster**)

swap labels and symbols when using **plot_lattice** to gain in homogeneity with **plot_base** and **plot_raster**

"histogram" legend were ok for **plot_lattice** but not for **plot_base** and **plot_raster** which were displaying `pch` instead of `lty`

'NA' is now an allowed value for `adjust_graph` argument in **data_rm** functions to force graph(s) removal

implement and export new **data_add_graphs** **data_rm_graphs** functions

**adjustGraphs** now internally uses graph as intput rather than indice(s) of graph(s) from `IFC_data` object (which is more versatile and allows for testing any graphs in the context of 'obj')

catch error in **adjustGraphs**

fix bugs when trying to adjust graph in several places

1 - when testing if drawable which was not working since about https://github.com/gitdemont/IFC/commit/030df2bfb25680e3159aa4a15bf0d550543d4188

2 - when applied GraphRegion(s) result(s) in not found pop(s)

3 - when trying to applyGatingStrategy which did not actually test if graph(s) could be produced

fix bug in **plotGraph**, **ExportToReport**, and **BatchReport** in feature summary computation for graphs containing overlay of populations (1D and 2D) for values of the region(s) drawn

- Fix **objectExtract**

fix the fix https://github.com/gitdemont/IFC/commit/34f33de1f916673855c2bc69c0c1e686bd608c58 which did not allow to correctly pass extra parameters to dots ...

- Fix **ExportToBatch**/**buildBatch**

absence of predefined report in the analysis file induced error when building the batch

fix bug when no compensation was provided

- FCS I/O

replace `first_only` argument in **readFCS** by `dataset` to allow user to choose which `dataset` to extract rather than only '1st' or 'All' (default behavior remains unchanged, however)

allow import/export of keywords

fix potential bug with extra keywords due to wrong regex identification

- Speed Gain

gray -> rgb color conversion improvement

use Rcpp::no_init when possible

modify internal **fastCbind** (the underlying hpp/cpp functions induced more overhead than R cbind)

- Misc

for progressbar, check if session missing and class rather than if length is 0

remove gc from fcs functions

allow to pass NULL and NA in **texttomatrix** / **addText** (before an error was thrown)

implement internal function **data_to_AST** to allow the creation of .ast file from `IFC_data` object

#### This leads to the following visible changes for the user
*See the bugs fix above*

*`first_only` argument is replaced by `dataset` for FCS reading*

*`IFC_data` object returned by ExtracFromFCS has gained an new `Keywords` member*

*new data_add_graphs and data_rm_graphs functions are now exported*

## 0.1.6
- CRAN submission

## 0.1.5.xxx
- File I/O

fix ExportToDAF that exports xml files containing NULL values

- issue with unexpected file names encoding

encode files names to native before passing them to underlying hpp/cpp files for reading

- object extraction

fix error in objectExtract when info argument is provided generated duplicated arguments

fix cpp_check_range that could lead to error when input can be summarize to a single value (this could happen when image object has no mask)

fix cpp_resize that did not take into account resizing information for masks

- XIF subset/merging

add `add_tracking` argument to mergeXIF, subsetXIF, XIFtoTIFF, and ExportToXIF functions to control whether exported file should contain tracking information

improve tracking of object_ID

improve XIF I/O speed (avoid conversion when building IFD and use raw instead thanks to new internal collapse_raw function)

- num to string

improve numeric to string conversion when using objects ids as names/dimnames (e.g. avoid 1e+05, use 100000 instead)

- population computation

create popsRename function to simplify population renaming in `IFC_data` object

fix bug with popsCompute for type "C" (combined) population computation when population definition contains repeated occurrence of same population name

#### This leads to the following visible changes for the user
*mergeXIF, subsetXIF and ExportToXIF functions gain a new `add_tracking` argument*

*a new popsRename function is now exported*


## 0.1.5
- CRAN release

fix trailing slash NOTE in last CRAN submission

fix cpp_normalize which did not take into account force_range parameter

## 0.1.4
- CRAN submission

## 0.1.2.101 to 0.1.3.104
- for CRAN

preparing for CRAN release

use doi instead of url

modify some links that have moved

- Improvements to graphs

improve drawing speed due to internal do.call that causes subtitute(x) to cast all x values to character

speed gain thanks to dissociation between `IFC_plot` object generation from plot and stats production

create internal plot_stats() dedicated function for stats extraction

create internal plot_lattice() dedicated function for plotting (also rename convert_to_baseplot to plot_base)

create internal plot_raster() and rasterplot() + cpp/hpp underlying functions to allow plot rasterization for very fast 2D graph drawing

better legend positioning

additional use of " " to prevent heightDetails potential error

fix bug with plotGraph which was opening a dev to retrieve lattice par settings

fix bug with grouping and precision = "light" when number of displayed population was <= 1

transformations (axes + density) should now be provided as character and are matched against allowed names

- Improvements to FCS file I/O

add text_only parameter

speed gain thanks creation of internal cpp/hpp functions after code profiling

better handle files with extra keywords segment

fix fcs import for multiple files merging when no spillover can be retrieved + when they have exactly same features names

fix large fcs file export

add support to bad-formatted fcs files (with empty data offsets in header)

warn user on missing additional datasets (i.e. non empty $NEXTDATA)

- Improvements to DAF / XIF file I/O

fix progress bar bug when extracting images/masks offsets

fix bug with population stats computation when IFC_data object contains tagged populations whose obj are not logical (i.e. numeric or integer)

fix bug when applying gating strategy with no graphs

better control/coercion on files read/write notably for graphs

- Misc

internally improve file scan by using raw vector target rather than string

fix internal whoami() which was erroring on namespaced calls (i.e. with :: or :::)

improve ExportToReport selection argument to allow graphs to be passed as a matrix that defines report layout

export new BatchReport function

#### This leads to the following visible changes for the user
*transformations should now be provided as character and will be matched against allowed names*

*a new BatchReport function is exported*

*readFCS gains a new 'text_only' entry in 'options' argument to allow user to only extract TXT segment*


## 0.1.2.100
- allow to pass `IFC_data` object in writeIFC()

- now handle `levelplot`

- improve .fcs input / output internally

speed gain and better accuracy notably for files stored with DATATYPE "I"

can handle files with non R standard bits length (i.e. not 1,2,4,8)

fix bug with large PnR e.g. "4294967296" which resulted in NA

- don't send error on partial mask retrieval

#### This leads to the following visible changes for the user
*buildGraph gains a new 'densitylevel' entry in BasePop controlling `levelplot` when 'type' is "density"*

*when fcs is imported @IFC_date, @IFC_file, @IFC_fileName, @IFC_dataset, @IFC_version, and @IFC_FCSversion are automatically filled. $FIL is not filled if found. $CYT is not filled if found except if user requires it.*

*readFCS gains a new 'force_header' entry in 'options' argument to allow user to control if data offset(s) should be taken from header only*


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
