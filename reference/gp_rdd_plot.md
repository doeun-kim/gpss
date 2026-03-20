# gp_rdd_plot

to draw an RD plot using the results obtained from gp_rdd()

## Usage

``` r
gp_rdd_plot(gp_rdd_res, l_col = "blue", r_col = "red")
```

## Arguments

- gp_rdd_res:

  a list-form results obtained from gp_rdd()

- l_col:

  a character value indicating the color of the left side of the cutoff
  point (default = "blue")

- r_col:

  a character value indicating the color of the right side of the cutoff
  point (default = "red")

## Value

A `ggplot` object showing the RD plot.

## Examples

``` r
library(ggplot2)
n <- 100
tau <- 3
cut <- 0
x <- rnorm(n, 0, 1)
y <- rnorm(n, 0, 1) + tau*ifelse(x>cut, 1, 0)
gp_rdd.out <- gp_rdd(x, y, cut)
gp_rdd_plot(gp_rdd.out) +
 geom_vline(xintercept = cut, lty="dashed")
```
