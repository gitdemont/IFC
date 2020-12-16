# NEWS

## 0.1.0
- fix bug when using verbose in Extract Imges/Masks + ExportTo/Display Gallery/Numpy family functions

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