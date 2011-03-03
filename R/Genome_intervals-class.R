# We define two classes for genome intervals. The basic class extends Interval_full with:
# a seq_name factor to represent chromosomes and more generally any sequence origin. The second also stores strand orientation.
# a logical vector "inter_base" (to flag zero-length features such as restriction sites) 


######## Class definitions
#### Genome_intervals

# extends Intervals_full of intervals package with a seqName factor

setClass(
		"Genome_intervals",
		representation = representation( annotation = "data.frame" ),
		prototype ( type ="Z"),
		contains = "Intervals_full",
		validity = function( object ) {
			# check 'annotation' data.frame
			if(!is.data.frame(object@annotation) || nrow( object@.Data ) != nrow( object@annotation ) )
				return("The 'annotation' slot must be a data.frame with as many rows as the endpoint matrix.")
			
			# Check 'seq_name' column within annotation data.frame slot
			if ( !('seq_name' %in% names(object@annotation)) || !is.factor( object@annotation$seq_name ) || any(is.na(object@annotation$seq_name)) )
				return( "The 'annotation' slot should have a column named 'seq_name' that is a factor and does not contain missing values." )
			
			# Check 'inter_base' column within annotation data.frame slot
			if ( !('inter_base' %in% names(object@annotation)) || !is.logical( object@annotation$inter_base ) || any(is.na(object@annotation$inter_base)) )
				return( "The 'annotation' slot should have a column named 'inter_base' that is logical and does not contain missing values." )
			
			if( object@type!='Z')
				return( "The intervals 'type' should be 'Z'." )
			
			return( TRUE )
		}
)

#### Genome_intervals_stranded

setClass(
		"Genome_intervals_stranded",
		contains = "Genome_intervals",
		validity = function( object ) {
			
			# Check 'strand' column within annotation data.frame slot
			if ( !('strand' %in% names(object@annotation)) || !is.factor( object@annotation$strand ) )
				return( "The 'annotation' slot should have column named 'strand' that is a factor." )
			
			if ( !nlevels(object@annotation$strand)==2 )
				return( "The 'strand' slot should be a factor with exactly two levels." )
			
			return( TRUE )
		}
)


######## Accessors and replacement methods


#### annotation
## implements the generics of package "BioBase"
setMethod(
		Biobase::annotation,
		signature( "Genome_intervals" ),
		function( object ) object@annotation
)
annotation <- Biobase::annotation

setMethod(
		Biobase::`annotation<-`, "Genome_intervals",
		function( object, value ) {
			object@annotation <- value
			return(object)
		}
)
`annotation<-` <- Biobase::`annotation<-`

#### inter_base
setGeneric( "inter_base", function(x) standardGeneric( "inter_base" ) )

setMethod(
		"inter_base",
		signature( "Genome_intervals" ),
		function( x ) x@annotation$inter_base
)

setGeneric( "inter_base<-", function( x, value ) standardGeneric( "inter_base<-" ) )

setReplaceMethod(
		"inter_base", "Genome_intervals",
		function( x, value ) {
			if(!( length( value ) %in% c(1,nrow(x)) ) | !is.logical(value) )
				stop( "The 'inter_base' argument should be a logical vector length equal to 1 or to the number of rows of the end point matrix." )
			if(length(value)==1) value = rep(value, nrow(x))
			x@annotation$inter_base <- value
			return(x)
		}
)


#### seq_name
setGeneric( "seq_name", function(x) standardGeneric( "seq_name" ) )

setMethod(
		"seq_name",
		signature( "Genome_intervals" ),
		function( x ) x@annotation$seq_name
)

setGeneric( "seq_name<-", function( x, value ) standardGeneric( "seq_name<-" ) )

setReplaceMethod(
		"seq_name", "Genome_intervals",
		function( x, value ) {
			if ( is.vector( value ) )
				value = factor(value)
			if(!( length( value ) %in% c(1,nrow(x)) ) )
				stop( "The 'seq_name' argument should be a vector or a factor of length equal to 1 or to the number of rows of the end point matrix." )
			if(length(value)==1) value = rep(value, nrow(x))
			x@annotation$seq_name <- value
			return(x)
		}
)

#### strand

setGeneric( "strand", function(x) standardGeneric( "strand" ) )

setMethod(
		"strand",
		signature( "Genome_intervals_stranded" ),
		function( x ) x@annotation$strand
)

setGeneric( "strand<-", function( x, value ) standardGeneric( "strand<-" ) )

setReplaceMethod(
		"strand", "Genome_intervals_stranded",
		function( x, value ) {
			if ( is.vector( value ) )
				value = factor(value)
			if(nlevels(value)!=2) 
				stop( "The 'strand' argument should be a vector with exactly 2 distinct values or a factor with 2 levels." )		
			if(!( length( value ) %in% c(1,nrow(x)) ) )
				stop( "The 'strand' argument should be a vector or a factor of length equal to 1 or to the number of rows of the end point matrix." )
			if(length(value)==1) value = rep(value, nrow(x))
			x@annotation$strand <- value
			return(x)
		}
)

#### type

setReplaceMethod(
		"type", "Genome_intervals",
		function( x, value ) {
			if ( length( value ) != 1 || value!=  "Z"  )
				stop( "The 'type' slot should be 'Z'." )
			x@type <- value
			return(x)
		}
)


######## Subsetting
setMethod(
		"[",
		signature( "Genome_intervals" ),
		function( x, i, j, ..., drop ) {
			if ( missing(i) ) i <- rep( TRUE, nrow(x) )
			if ( missing(j) ) {
				# Preserve class. Note that both [i,] and [i] syntax subset rows.
				if ( is.character(i) ) i <- match( i, rownames( x ) )
				x@annotation <- x@annotation[i,,drop=FALSE]
			}
			callNextMethod( x, i, j, ..., drop )
		}
)


setMethod(
		"[<-",
		signature( x = "Genome_intervals", i = "ANY", j = "missing", value = "Genome_intervals" ),
		function( x, i, j, value ) {
			#### Error checking
			if ( is.character(i) ) i <- match( i, rownames( x ) )
			if ( any( is.na( i ) ) )
				stop( "Cannot assign to NA indices or row names which do not exist." )
			n <- length( (1:nrow(x))[i] )
			if ( n != nrow( value ) )
				stop( "Replacement object is of the wrong size." )            
			
			if( length(annotation(value)) != length(annotation(x)) )
				stop("Number of columns of annotation slots do not match. Check if you're assigning a Genome_intervals_stranded into rows of a Genome_intervals or vice-versa.")
			if( any( names(annotation(value)) != names(annotation(x)) ) )
				stop("Names of annotation do not match. Check if you're assigning a Genome_intervals_stranded into rows of a Genome_intervals or vice-versa.")
			
			#### Intervals
			x@.Data[i,] <- value@.Data
			x@closed[i,] <- value@closed
			
			## Annotation            
			annotation(x)[i,] <- annotation(value)
			
			#### Rownames
			has_names_x <- !is.null( rownames(x) )
			has_names_value <- !is.null( rownames(value) )
			if ( has_names_x ) {
				if ( has_names_value ) rownames(x)[i] <- rownames(value)
				else rownames(x)[i] <- ""
			}
			else {
				if ( has_names_value ) {
					rownames(x) <- rep( "", nrow(x) )
					rownames(x)[i] <- rownames(value)
				}
			}
			return(x)
		}
)


setMethod("$", "Genome_intervals", function(x, name) {
			eval(substitute(annotation(x)$NAME_ARG, list(NAME_ARG=name)))
		})

setReplaceMethod("$", "Genome_intervals", function(x, name, value) {
			x[[name]] <- value
			x
		})

setMethod("[[", "Genome_intervals", function(x, i, j, ...) annotation(x)[[i]] )

setReplaceMethod("[[",
		signature=signature(x="Genome_intervals"),
		function(x, i, j, ..., value) {
			annotation(x)[[i]] <- value
			x
		})

######## Coercion

setMethod(
		"coerce",
		signature( from = "Genome_intervals", to = "Intervals_full" ),
		function( from, to, strict ) {
			new(
					"Intervals_full",
					from@.Data,
					type = type( from ),
					closed = closed(from)
			)
		}
)

setMethod(
		"coerce",
		signature( from = "Genome_intervals", to = "character" ),
		function( from, to, strict ) {
			if ( nrow( from ) == 0 )
				return( character() )
			else {
				# call to Intervals coercion method for the intervals
				ints <- as( as(from, "Intervals_full"), "character")
				# add seq_name in first column, inter_base in last column
				result <- paste(seq_name(from), ints, ifelse( inter_base(from), "*", ""))
				return( result )
			}
		}
)

setMethod(
		"coerce",
		signature( from = "Genome_intervals_stranded", to = "character" ),
		function( from, to, strict ) {
			if ( nrow( from ) == 0 )
				return( character() )
			else {
				# call to Intervals coercion method for the intervals
				ints <- as( as(from, "Intervals_full"), "character")
				# add seq_name and strand in first columns, inter_base in last column
				result <- paste(seq_name(from), strand(from), ints, ifelse( inter_base(from), "*", "") )
				return( result )
			}
		}
)


## more user-friendly constructor function:
GenomeIntervals <- function(chromosome, start, end, strand=NULL,
                            inter.base=NULL, leftOpen=NULL,
                            rightOpen=NULL, ...)
{
  stopifnot(length(chromosome)==length(start),
            length(chromosome)==length(end))
  ## any inter-base intervals specified
  if (is.null(inter.base)){
    inter.base <- vector("logical", length(chromosome))
  } else {
    stopifnot(is.logical(inter.base), length(inter.base)==length(chromosome))
  }
  ## by default, all intervals are assumed to be closed,
  ##  unless specified otherwise:
  if (!is.null(leftOpen))
    stopifnot(is.logical(leftOpen), length(leftOpen)==length(chromosome))
  else
    leftOpen <- vector("logical", length(chromosome))
 if (!is.null(rightOpen))
    stopifnot(is.logical(rightOpen), length(rightOpen)==length(chromosome))
  else
    rightOpen <- vector("logical", length(chromosome))
  
  ## prepare annotation data.frame for object:
  annoDf <- data.frame("seq_name"=factor(chromosome),
                       "inter_base"=inter.base)
  ## additional annotation columns?
  further.args <- as.list(match.call(expand.dots=FALSE)[["..."]])
  if (length(further.args)>0){
    for (i in 1:length(further.args))
      annoDf[[names(further.args)[i]]] <- eval(further.args[[i]])
  }

  ## create object depending on whether or not the strand was specified
  if (!is.null(strand)){
    ## create object of class "Genome_intervals_stranded"
    stopifnot(length(strand)==length(chromosome),
              all(strand %in% c("+", "-")))
    annoDf$strand <- factor(strand, levels=c("+", "-"))
    obj <- new("Genome_intervals_stranded",
               cbind(start, end),
               closed=cbind(!leftOpen, !rightOpen),
               annotation=annoDf)
  } else {
    ### non-stranded intervals
    obj <- new("Genome_intervals",
               cbind(start, end),
               closed=cbind(!leftOpen, !rightOpen),
               annotation=annoDf)
  }
  stopifnot(validObject(obj))
  return(obj)
}# GenomeIntervals

### sort method for Genome_intervals objects:
setMethod("sort", signature(x="Genome_intervals"),
          function(x, decreasing=FALSE, ...) {
            chr <- gsub("^[Cc]hr", "",
                        as.character(seq_name(x)))
            ## replace some chromosome names by numbers for sorting:
            chr <- gsub("X$","100", chr)
            chr <- gsub("Y$","200", chr)
            chr <- gsub("MT?$","300", chr)
            chr <- gsub("_random$","000", chr)
            suppressWarnings(chr <- as.numeric(chr))
            ## are there still non-numeric entries left in chr?
            if (any(is.na(chr))) chr <- as.character(seq_name(x))
            ord <- order(chr, x[,1], x[,2], decreasing=decreasing)
            return(x[ord])
} ) # sort
