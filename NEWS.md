# n1qn1 6.0.1-12.9000 (development version)

* Verify no upstream Scilab Fortran changes need porting: the two 2020 Scilab
  commits (`977d4b6c`, `b9c2de06`) fixed Fortran 77 rank-mismatch errors in
  `n1qn1a.f` for gfortran 10 compatibility by declaring `f` as `f(1)`. The
  f2c-translated C code already uses a separate scalar `fb` variable, so these
  changes are already correctly represented in the C translation.

* Make package thread-safe: convert global state to `thread_local`, remove
  `static` from local variables in Fortran-translated C code, add integer
  overflow guards, and fix memory leaks in the R callback wrappers.

* Review all C integer types; confirmed existing types are correct given
  Fortran ABI constraints (`int*` required for Fortran INTEGER) and R's
  `INTEGER()` returning `int*`.

* Fix undefined behavior: add virtual destructor to `EvalBase` so deleting
  derived callback objects through a base pointer is well-defined.

* Fix memory leak: validate `n` and `nzm` before allocating callback objects,
  preventing a leak when `Rf_error` is called on bad inputs.

* Fix `restart=TRUE`: typo caused mode to remain 2 instead of being set to 3.

* Fix `assign=TRUE`: wrong field name meant the compressed Hessian was never
  stored in the supplied environment.

* Add tests for restart, assign, single-variable optimization, repeated calls,
  and input error conditions.

* Make package thread-safe: convert global state to `thread_local`, remove
  `static` from local variables in Fortran-translated C code, add integer
  overflow guards, and fix memory leaks in the R callback wrappers.

* Review all C integer types; confirmed existing types are correct given
  Fortran ABI constraints (`int*` required for Fortran INTEGER) and R's
  `INTEGER()` returning `int*`.

# n1qn1 6.0.1-12

* Add non binary (function pointer) interface/api so that `nlmixr2est`
  will not have to be re submitted when changes to this package occur.

# n1qn1 6.0.1-11

* Added strict prototype fixes as requested by CRAN

* Added a `NEWS.md` file to track changes to the package.
