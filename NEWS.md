## News for Package `DMMF`

#### Changes in DMMF version 0.5.2.0
- Remove dependencies of retiring packages (rgdal, rgeos, maptools) 

#### Changes in DMMF version 0.5.1.2
- Bug fix for Sinkfill algorithm

#### Changes in DMMF version 0.5.1.1
- Adopt (x=0.0d0 and NaN = 0.0d0/x) for NaN value

#### Changes in DMMF version 0.5.1.0
- Adopt `ieee_arithmetic` for NaN value

#### Changes in DMMF version 0.5.0.4
- Change NaN value from NaN to -999999 in Fortran code.

#### Changes in DMMF version 0.5.0.3
- Minor update to debug compile errors by modifying `DMMF.h` and `R_init_DMMF.c` (remove lines related to all Fortran functions in both files).

#### Changes in DMMF version 0.5.0.2
- Modifying SinkFill algorithm to give correct results.

#### Changes in DMMF version 0.5.0.1
- From the bug reporting from Mohammad Fereshtehpour, we modified codes to compare xres and yres of DEM from "xres(DEM) == yres(DEM)" to "all.equal(xres(DEM), yres(DEM))".

#### Changes in DMMF version 0.4.5.1
- Add function `DMMFc` which is useful when input maps are categorical.

#### Changes in DMMF version 0.4.0.1
- Modifying codes using `isnan` from GCC function in accordance with fortran 95 standard.
- Change in raster maps in `Potato.Convex` and `Potato.Concave` to have inMemory options.

#### Changes in DMMF version 0.3.1.0
- Modifying function `SinkFill` to create the sinkfilled map and the partition of the catchment. 

#### Changes in DMMF version 0.3.0.0
- Add function `MapChecker` for detect boundary, sink, and stand cells.
- `B_checker` function is merged to the `MapChecker` function.
- Modify the `SinkFill` algorithm and `BoundaryCheck` algorithm to work properly.

#### Changes in DMMF version 0.2.8.1
- Add `B_checker` function to detect boundary cells of the DEM and manual.
- Modify the `SinkFill` algorithm and `BoundaryCheck` algorithm to work properly.
- It is not certain but there is a problem when using the values of NA from R in Fortran.

#### Changes in DMMF version 0.2.8.0
- Resolve the note message when build the package in the windows systems
- Error message is that `"Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols"`.
- Follow the `init.c` and header files from `stats` package included in R.
- Add prefix "C\_" for Fortran functions when they are called in R as `.Fortran()`.
- Ready to submit package to CRAN.

#### Changes in DMMF version 0.2.7.4
- Modify `/src/init.c` for routines to have proper types of elements.

#### Changes in DMMF version 0.2.7.1
- Add `/src/init.c` to avoid the error message under developement version using `tools::package_native_routine_registration_skeleton('DMMF')`.
- Before using `package_native_routine_registration_skeleton`, working directory should be set as parent directory of the package.

#### Changes in DMMF version 0.2.7
- Change `NAMESPACE` to avoid note messages when build in unstable developement version of R in win-builder.
- Remove line `exportPattern("^[[:alpha:]]+")`
- Add line `export("DMMF", "DMMF_Simple", "SinkFill")`

#### Changes in DMMF version 0.2.5
- Change `NAMESPACE` from `useDynLib(DMMF)` to `useDynLib(DMMF, .registration = TRUE)` to remove note messages when build in unstable developement version of R in win-builder.

#### Changes in DMMF version 0.2.4
- Add `\dontrun{}` option for `SinkFill` and `DMMF_simple` example to avoid error messages in win-builder. 

#### Changes in DMMF version 0.2.4
- Add `\dontrun{}` option for `DMMF` example to avoid error messages in win-builder. 

#### Changes in DMMF version 0.2.3
- Modifying the source code of MMF Algorithm from `SW_t = dmax1( 0.0, SW - SW_fc )` to `SW_t = dmax1( 0.0d0, SW - SW_fc )` to make 0.0 as double.

#### Changes in DMMF version 0.2.2
- Modifying NAMESPACE as import `rgdal` package (import("rgdal"))and import `is` from `methods` package (importfrom("methods", "is")).

#### Changes in DMMF version 0.2.1
- Bug fixed: All non-ASCII characters in the manuals are corrected to ASCII character to remove warnings when check the package with `R CMD CHECK --as-CRAN`.
- Bug fixed: Title and package description parts of the `/man/DMMF-package.Rd` is changed to remove warnings when check the package with `R CMD CHECK --as-CRAN`.
- README.md and NEWS.md files are added to the package.

#### Changes in DMMF version 0.2.0
- The new function named `DMMF_Simple` is added to the package which is the simpler version of `DMMF` function.
- The `/src/2_2_BoundaryChecker.f95` is changed to convert logical matrix to real value matrix of 1 and 0 smoothly without warning message of "implicit conversion from logical to integer".

#### Changes in DMMF version 0.1.0
- The manual of functions in the package are added to the package.
- Datasets from two potato fields (`Potato.Concave` and `Potato.Convex`) are added to the package.


