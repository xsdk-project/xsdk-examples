# PLASMA examples
Example codes demonstrating the use of PLASMA using other XSDK packages.

## Using SLATE and its compute layers
`ex1solve` solves a system of linear equations by using Level 3 BLAS from
SLATE's BLAS++ instead of PLASMA's internal interface. The example accepts an
optional command line  parameter to specify linear system size like so:
```
./ex1solve [--n=1000]
```
If the size is not specified then 1000 is used.
