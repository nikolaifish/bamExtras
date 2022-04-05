# bamExtras
`bamExtras` is an R package to that contains functions and data to support
stock assessments of fish stocks in the Southeast US Atlantic using the Beaufort
Assessment Model (BAM)

# Example usage

## Load the package

First, install and load the package. 

```{r, echo=TRUE, message=FALSE}
# Install and load package
devtools::install_github("nikolaifish/bamExtras", build_vignettes = TRUE)
library(bamExtras)
```


## Vignettes available
Examples of how this package can be used are provided in the vignettes (i.e. long-form guides).

The `intro` vignette introduces the basics of the package.
```R
vignette("intro","bamExtras")
```
