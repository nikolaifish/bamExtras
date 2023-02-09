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

# Disclaimer
“This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.”
