#!/usr/bin/env Rscript
library(testthat)
library(SWATH2stats)
all <- try(testthat::test_dir("tests/testthat", reporter="summary"))
summary(all)
