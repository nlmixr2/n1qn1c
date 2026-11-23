# This gives the function pointers in the n1qn1 library

Using this will allow C-level linking by function pointers instead of
abi.

## Usage

``` r
.n1qn1ptr()
```

## Value

list of pointers to the n1qn1 functions

## Author

Matthew L. Fidler

## Examples

``` r

.n1qn1ptr()
#> $n1qn1F
#> <pointer: 0x70b0d2c57990>
#> 
#> $n1qn1F2
#> <pointer: 0x70b0d2c57a00>
#> 
#> $n1qn1_
#> <pointer: 0x70b0d2c597d0>
#> 
```
