library(knitr)
library(gRbase)
library(Rgraphviz)
library(gRim)
library(unifips)
library(ggplot2)
library(rbenchmark)
library(parallel)
library(Matrix)

devtools::load_all("unifips")
source("script/multiplot.R")
source("script/bench-fun.R")

make_initc("unifips")
compileAttributes("unifips")
devtools::load_all("unifips")
ggmu(~.^., data=mtcars)

