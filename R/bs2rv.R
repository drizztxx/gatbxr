#' @title Binary string to real vector
#'
#' @description
#' This function decodes binary chromosomes into vectors of reals. The
#' chromosomes are seen as the concatenation of binary strings of given
#' length, and decoded into real numbers in a specified interval using
#' either standard binary or Gray decoding.
#'
#' @param Chrom matrix containing the chromosomes of the current
#' population. Each row corresponds to one
#' individual's concatenated binary string
#' representation. Leftmost bits are the most significant bits and
#' rightmost are the least significant bits.
#' @param FieldD a list describing how to decode each substring in the
#' chromosome. It must contain the following 7 elements:
#' \describe{
#'   \item{\code{prec}}{a number giving the length (precision) of each substring;}
#'   \item{\code{lb}, \code{ub}}{lower and upper bounds of each variable;}
#'   \item{\code{code}}{\code{"binary"} or \code{"gray"}, indicating how each
#'   substring is to be decoded;}
#'   \item{\code{scale}}{\code{"arith"} or \code{"log"}, indicating whether to use
#'   arithmetic or logarithmic scaling;}
#'   \item{\code{lbin}, \code{ubin}}{logical values indicating whether or not to
#'   include each bound in the representation range.}
#' }
#' Except for \code{prec}, each element may be a single value applied to all
#' variables or a vector with one value per variable.
#' @return a real matrix containing the population phenotypes. Each row
#' corresponds to one individual and each column to one variable.
#' @export
#' @author
#' The original 'MATLAB' implementation of bs2rv was written by Carlos Fonseca and
#' updated by Andrew Chipperfield.
#' The R implementation was written by Liang Zhao.
#' @seealso \code{\link{crtbp}}
#' @example examples/bs2rv.R

bs2rv <- function(Chrom, FieldD){
  ## Identify the population size (Nind) and the chromosome length (Lind)
  Chrom <- as.matrix(Chrom)
  Nind <- NROW(Chrom); Lind <- NCOL(Chrom)

  ## Check FieldD
  fields <- c("prec", "lb", "ub", "code", "scale", "lbin", "ubin")
  if (!is.list(FieldD) || !all(fields %in% names(FieldD))) {
    stop("FieldD must be a list containing the elements ",
         paste(fields, collapse = ", "))
  }

  ## Identify the number of decision variables (Nvar)
  prec <- FieldD$prec
  if (length(prec) != 1) stop("prec must be a scalar")
  if (Lind %% prec) stop("prec disagrees with chromosome length")
  Nvar <- Lind/prec

  ## Get substring properties, one value per variable
  expand <- function(x, name){
    if (length(x) == 1) x <- rep(x, Nvar)
    if (length(x) != Nvar) stop(name, " must have length 1 or ", Nvar)
    x
  }
  lb    <- expand(FieldD$lb, "lb")
  ub    <- expand(FieldD$ub, "ub")
  code  <- expand(FieldD$code, "code")
  scale <- expand(FieldD$scale, "scale")
  lin   <- expand(as.logical(FieldD$lbin), "lbin")
  uin   <- expand(as.logical(FieldD$ubin), "ubin")

  ## Check substring properties for consistency
  if (!all(code %in% c("binary", "gray"))) stop("code must be \"binary\" or \"gray\"")
  if (!all(scale %in% c("arith", "log"))) stop("scale must be \"arith\" or \"log\"")
  islog <- scale == "log"
  if (any(lb[islog] * ub[islog] <= 0)) {
    stop("Log-scaled variables must not include 0 in their range")
  }

  ## Decode chromosomes
  logsgn <- sign(lb)
  lb[islog] <- log(abs(lb[islog]))
  ub[islog] <- log(abs(ub[islog]))
  delta <- ub - lb

  Prec <- 0.5^prec
  num <- (!lin) * Prec
  den <- (lin + uin - 1) * Prec

  Phen <- matrix(0, Nind, Nvar)
  for (i in seq_len(Nvar)){
    idx <- ((i-1)*prec+1):(i*prec)
    Sub <- Chrom[, idx, drop = FALSE]
    if (code[i] == "gray" && prec > 1){ ## Gray decoding
      Sub <- t(apply(Sub, 1, cumsum)) %% 2
    }
    Phen[, i] <- Sub %*% (0.5^(1:prec))
    Phen[, i] <- lb[i] + delta[i] * (Phen[, i] + num[i]) / (1 - den[i])
  }
  if (any(islog)){
    Phen[, islog] <- rep(logsgn[islog], each = Nind) * exp(Phen[, islog])
  }
  return(Phen)
}
