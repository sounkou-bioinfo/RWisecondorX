# Extract fractions suitable for regression

For `CombinedControlGroup`: 22 x N matrix. For `SeparatedControlGroup`:
44 x N matrix (rows "1F".."22F","1R".."22R").

## Usage

``` r
fractions_for_regression(cg, ...)
```

## Arguments

- cg:

  A `NIPTControlGroup` object.

- ...:

  Reserved for S7 method dispatch; currently unused.

## Value

A numeric matrix.
