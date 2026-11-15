## Find the function given by the operator argument (REC_F, SEL_F or MUT_F) of
## a high-level function. The operators of this package are found even when it
## is not attached (e.g. when called as gatbxr::recombin() from another
## package); any other name is looked up from 'envir', the environment the
## high-level function was called from, so user-defined operators work too.
get_operator <- function(FUN, arg, envir){
  if (is.function(FUN)) return(FUN)
  if (!is.character(FUN) || length(FUN) != 1 || is.na(FUN)) {
    stop(arg, " must be a string of function name", call. = FALSE)
  }
  fun <- get0(FUN, envir = asNamespace("gatbxr"), mode = "function", inherits = FALSE)
  if (is.null(fun)) fun <- get0(FUN, envir = envir, mode = "function")
  if (is.null(fun)) stop(arg, ": could not find function \"", FUN, "\"", call. = FALSE)
  fun
}
