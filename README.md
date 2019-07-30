---
title: "gglogo"
author: "Eric Hare, Heike Hofmann"
date: "July 30, 2019"
output: 
  html_document:
    keep_md: true
---



R package for creating sequence logo plots

[![CRAN Status](http://www.r-pkg.org/badges/version/gglogo)](https://cran.r-project.org/package=gglogo) [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/gglogo)](http://www.r-pkg.org/pkg/gglogo) 
[![Travis-CI Build Status](https://travis-ci.org/heike/gglogo.svg?branch=master)](https://travis-ci.org/heike/gglogo)

# Installation

`gglogo` is available from CRAN (version 0.0.2):

```r
install.packages("gglogo")
```


The development version is available from Github:


```r
# install.packages("devtools")
devtools::install_github("heike/gglogo", build_vignettes = TRUE)
```

# Getting Started

Load the library


```r
library(gglogo)
```

Load a dataset


```r
data(sequences)
```

## A first sequence logo plot


```r
library(ggplot2)

ggplot(data = ggfortify(sequences, peptide)) +      
  geom_logo(aes(x = position, y = bits, group = element, 
     label = element, fill = interaction(Polarity, Water)),
     alpha = 0.6)  +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "bottom")
```

![](README_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

