#' renewr: Recovery of modulated renewal processes via Gaussian process modelling
#'
#' @description A compact toolkit for Gaussian-process modelling of time-series discrete event data
#' \if{html}{\figure{hex_sticker_RenewR.png}{options: style="float: right" width="150" alt="Hex Sticker"}}
#' @author Caius Gibeily \email{cgibeil@emory.edu}
#' @name renewr-package
#' @aliases renewr
#' @useDynLib renewr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import rkriging
#' @import tibble
#' @import flexsurv
#' @import rTensor
#' @import ggplot2
#' @import ggdist
#' @import patchwork
#' @import loo
#' @importFrom rstan sampling rstan_options extract
#' @importFrom data.table CJ
#' @importFrom cowplot plot_grid
#' @importFrom RColorBrewer brewer.pal
#' @importFrom MASS mvrnorm
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.7. https://mc-stan.org
"_PACKAGE"
NULL
