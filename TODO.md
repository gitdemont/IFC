# TODO
- add object tracking to ExportToXIF, include origin (fileName + checksum + object_id)
- use parallelization a good candidate would be cpp_extract (but how to pass multiple arguments in a threadsafe manner to std::transform in RcppParallel::parallelFor ?)
- add getAborted() to easily retrieve files for wich batch has failed
- extract software + version information in getInfo(), ExtractFromDAF(), ExtracFromXIF()
