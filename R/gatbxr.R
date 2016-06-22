#' @title Gentic Algorithm Toolbox Implemented by R
#' @description
#' R version of genetic toolbox for implementing a wide range of genetic algorithm methods.
#' @details
#' To create populations:
#' \itemize{
#'  \item{\code{\link{crtbase}}} Create a base vector
#'  \item{\code{\link{crtbp}}} Create arbitrary discrete random populations
#'  \item{\code{\link{crtrp}}} Create real-valued initial population
#' }
#' 
#' Fitness assignment:
#' \itemize{
#'  \item{\code{\link{ranking}}} Generalised rank-based fitness assignment
#'  \item{\code{\link{scaling}}} Proportional fitness scailing
#' }
#' 
#' Selection functions:
#' \itemize{
#'  \item{\code{\link{reins}}} Uniform random and fitness-based reinsertion
#'  \item{\code{\link{rws}}} Roulette wheel selection
#'  \item{\code{\link{select}}} High-level selection routine
#'  \item{\code{\link{sus}}} Stochastic universal sampling
#' }
#' 
#' Mutation operators:
#' \itemize{
#'  \item{\code{\link{mut}}} Discrete mutation
#'  \item{\code{\link{mutate}}} High-level mutation function
#'  \item{\code{\link{mutbga}}} Real-value mutation
#' }
#' 
#' Crossover operators:
#' \itemize{
#'  \item{\code{\link{recdis}}} Discrete recombination
#'  \item{\code{\link{recint}}} Intermediate recombination
#'  \item{\code{\link{reclin}}} Line recombination
#'  \item{\code{\link{recmut}}} Line recombination with mutation features
#'  \item{\code{\link{recombin}}} High-level recombination operator
#'  \item{\code{\link{xovdp}}} Double-point crossover
#'  \item{\code{\link{xovdprs}}} Double-point reduced surrogate crossover
#'  \item{\code{\link{xovmp}}} General multi-point crossover
#'  \item{\code{\link{xovsh}}} Shuffle crossover
#'  \item{\code{\link{xovshrs}}} Shuffle reduced surrogate crossover
#'  \item{\code{\link{xovsp}}} Single-point crossover
#'  \item{\code{\link{xovsprs}}} Single-point reduced surrogate crossover
#' }
#' 
#' Subpopulation support:
#' \itemize{
#'  \item{\code{\link{migrate}}} exchange individuals between subpopulations
#' }
#' 
#' Utility functions:
#' \itemize{
#'  \item{\code{\link{bs2rv}}} binary string to real-value conversion
#' }
#' @examples
#' ## Simple Genetic Algorithm
#' @example examples/sga.R
#' @name gatbxr-package
#' @aliases gatbxr
#' @author David Zhao \email{wethenwethen@gmail.com}
#' @docType package
NULL