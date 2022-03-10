#Thanks to Tycho Bismeijer for implementation
requireNamespace('rmarkdown')

p <- local({
  p <- list()
  args <- commandArgs(T)
  stopifnot(length(args) >= 2)
  p$rmd <- args[1]
  p$out  <- args[2]
  
  if (length(args) > 2) {
    params <- args[3:length(args)]
    stopifnot(length(params) %% 2 == 0)
    params_names <- params[seq(1, length(params), 2)]
    params_values <- params[seq(1, length(params), 2)+1]
    stopifnot(stringr::str_sub(params_names, 1, 2) == '--')
    p$params <- structure(as.list(params_values),
                          names=stringr::str_sub(params_names, 3, -1))
  } else {
    p$params=NULL
  }
  p
})

#rmarkdown::render(p$rmd, rmarkdown::md_document(variant='markdown'), p$out, params=p$params)
rmarkdown::render(input = p$rmd, output_file = p$out, params=p$params)
