##' c extension for the \code{genomeIntervals package}
##'
##' This function combines severl genome intervals (stranded or not)
##' objects into a single one.
##'
##' If the arguments have mixed classes ( both \code{\linkS4class{Genome_intervals}}
##' or \code{\linkS4class{Genome_intervals_stranded}}), then
##' they are coerced to \code{\linkS4class{Genome_intervals}} before combination.
##' Otherwise, the common class is used. Ifa list is provided with \code{NULL}
##' entries, these are discarded. If a vector of object is provided with non \code{genomeIntervals}
##' classes, then a list, ordered as the input vector, is returned.
##'
##' @aliases c c,Genome_intervals-methods
##' c,Genome_intervals_stranded-methods
##' @name c
##' @rdname c.Genome_intervals
##' @param ... two or more \code{\linkS4class{Genome_intervals}} or \code{\linkS4class{Genome_intervals_stranded}} objects.
##' @param recursive inherited from the base \{c} function definition and not used.
##' @return
##' \itemize{
##'   \item A single \code{\linkS4class{Genome_intervals}} or
##'   \code{\linkS4class{Genome_intervals_stranded}} object. Input objects are
##'   combined in their order of appearance in the the argument list.
##'   \item If any input argument is not a \code{\linkS4class{Genome_intervals}}, \code{list(...)} is
##'   returned instead.
##' }##'
##' @examples
##'   ##' load toy examples
##'   data("gen_ints")
##'
##'   ##' combine i and j returns a Genome_intervals_stranded object
##'   c( i, j )
##'
##'   ##' combine a not-stranded and a stranded returns a not-stranded object
##'   c( as(i, "Genome_intervals"), j )
##'
setMethod("c",
          signature=("Genome_intervals"),
          definition=function(...,recursive=FALSE
          ){
            args <- list(...)

            ## Drop NULL arguments
            if ( any( sapply( args, is.null ) ) )
              args <- args[ !sapply( args, is.null ) ]

            ## If anything else than a Genome_intervals (stranded) returns a list
            classes <- sapply( args, class )
            if ( !all( classes %in% c( "Genome_intervals", "Genome_intervals_stranded" ) ) )
              return( list( ... ) )

            ## If mixed Genome_intervals types (stranded and not), coerce to not stranded
            same_class <- all( classes == "Genome_intervals" )
            if ( !same_class  ) {
              warning( "Coercion to 'Genome_intervals' required.", call. = FALSE )
              return( do.call( c, lapply( args, as, "Genome_intervals" ) ) )
            }

            # rbinds the data frame if possible
            annot =  try( do.call( rbind, lapply( args, function(x) annotation(x) ) ), silent = TRUE )
            if( class(annot) == "try-error")
              stop("impossible to combine the annotation slots:\n", annot, call.=FALSE)

            # actual appending
            result <- new(
              "Genome_intervals",
              do.call( c, lapply( args, as, "Intervals_full" ) ),
              annotation = annot
            )
            return( result )
          }
)

setMethod("c",
          signature=("Genome_intervals_stranded"),
          definition=function(...,recursive=FALSE
          ){
            args <- list(...)

            ## Drop NULL arguments
            if ( any( sapply( args, is.null ) ) )
              args <- args[ !sapply( args, is.null ) ]

            ## If anything else than a Genome_intervals (stranded) returns a list
            classes <- sapply( args, class )
            if ( !all( classes %in% c( "Genome_intervals", "Genome_intervals_stranded" ) ) )
              return( list( ... ) )

            ## If mixed Genome_intervals types (stranded and not), coerce to not stranded
            same_class <- all( classes == "Genome_intervals_stranded" )
            if ( !same_class  ) {
              warning( "Coercion to 'Genome_intervals' required.", call. = FALSE )
              return( do.call( c, lapply( args, as, "Genome_intervals" ) ) )
            }

            # rbinds the data frame if possible
            annot =  try( do.call( rbind, lapply( args, function(x) annotation(x) ) ), silent = TRUE )
            if( class(annot) == "try-error")
              stop("impossible to combine the annotation slots:\n", annot)

            # actual appending
            result <- new(
              "Genome_intervals_stranded",
              do.call( c, lapply( args, as, "Intervals_full" ) ),
              annotation = annot
            )
            return( result )
          }
)