
#' Wrapper for Rmarkdown parametrization
#' 

library(methods)
library(stringr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

args_array <- unlist(strsplit(args, " "))

capture_args <- function(args_array) {
  name_pos <- which(stringr::str_detect(args_array, "^\\-\\-"))
  arg_pos <- setdiff(name_pos + 1, name_pos)
  arg_pos <- arg_pos[arg_pos <= length(args_array)]
  
  args <- as.list(args_array[arg_pos])
  arg_names <- stringr::str_replace(args_array[name_pos], "^\\-\\-", "")
  
  stopifnot(!any(duplicated(arg_names)))
  
  names(args) <- arg_names
  return(args)
}

arg_list <- capture_args(args_array)

stopifnot(all(c('input_file') %in% names(arg_list)))

input_rmd_file <- arg_list[['input_file']]

other_param_names <- setdiff(names(arg_list), c("input_file"))
other_params <- arg_list[other_param_names]
other_param_string <- paste(lapply(names(other_params), function(x) {
  paste0(x, "=", '"', other_params[[x]], '"')
}), collapse = ", ")

## Test existence of rmarkdown file
stopifnot(file.exists(input_rmd_file))

cmd <- sprintf('cp %s %s/; Rscript -e \'rmarkdown::render("%s", "html_notebook", params = list(%s))\'', input_rmd_file, getwd(), paste0(getwd(), "/", basename(input_rmd_file)), other_param_string)
print(cmd)

system(cmd)

cat("Completed.\n")


