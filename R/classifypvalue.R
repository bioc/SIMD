#' calculate P-value in code EMtest.
#'
#' @param type1 The first colum of the first matrix.
#' @param type2 The second colum of the first matrix.
#' @param type3 The third colum of the first matrix.
#' @param type4 The fourth colum of the first matrix.
#' @param sm1chring1 The first colum of the second matrix.
#' @param sm1chring2 The second colum of the second matrix.
#' @param sm1chring3 The third colum of the second matrix.
#' @param sm1chring4 The forth colum of the second matrix.
#' @param p P-value.
#' @param typelength The nrows of the first matrix.
#' @param sm1chringlength The nrows of the second matrix.
#' @param pvalue A vector, the length equals to the nrows of the second matrix.
#' @return The probability.
#' @keywords internal
#' @useDynLib SIMD pvalueclassify

classifypvalue <- function(type1, type2, type3, type4, 
                           sm1chring1, sm1chring2, sm1chring3, 
                           sm1chring4, p, typelength, sm1chringlength,
                           pvalue=rep(0,length(sm1chring1))){
    problity <- .C("pvalueclassify", as.integer(type1), 
                   as.integer(type2), as.integer(type3), 
                   as.integer(type4), as.integer(sm1chring1),
                   as.integer(sm1chring2), as.integer(sm1chring3), 
                   as.integer(sm1chring4), as.double(p), 
                   as.integer(typelength), as.integer(sm1chringlength),
                   as.double(pvalue))
    return(problity[[12]])
}

.onUnload <- function (libpath) {library.dynam.unload("SIMD", libpath)}







