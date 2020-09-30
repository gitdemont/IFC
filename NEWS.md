# NEWS

## 0.1.0
- fix bug when using verbose in Extract Imges/Masks + ExportTo/Display Gallery/Numpy family functions

- fix bug in objectDisplay when image has more then one class

- fix bug in ExportToNumpyobject size was not correctly set

- pass all cpp functions to header to allow more portability (side effect package installation is faster)

- reintegrate num_to_string function using format() did not produce desired results

## 0.0.9
- allow graph conversion from lattice to base

- modify plotGraph return object to allow to pass it to base plot + better handle graphics device

- fix solaris installation error