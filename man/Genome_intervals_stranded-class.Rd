\name{Genome_intervals_stranded-class}
\docType{class}

\alias{Genome_intervals_stranded-class}
\alias{coerce,Genome_intervals_stranded,character-method}

\alias{strand}
\alias{strand<-}
\alias{strand,Genome_intervals_stranded-method}
\alias{strand<-,Genome_intervals_stranded-method}

\title{Class "Genome\_intervals\_stranded" }
\description{A set of genomic intervals with a specified strand.}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{See \code{\linkS4class{Genome_intervals} }}
    \item{\code{annotation}:}{A \code{data.frame} (see \code{\linkS4class{Genome_intervals}} for basic requirements).
    The \code{annotation} moreover has a \code{strand} column that is a factor with exactly two levels(typically \code{"+"} and \code{"-"}). }
    \item{\code{closed}:}{See \code{\linkS4class{Genome_intervals} } }
    \item{\code{type}:}{See \code{\linkS4class{Genome_intervals} } }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Genome_intervals}"}, directly.
Class \code{"\linkS4class{Intervals_full}"}, by class "Genome\_intervals", distance 2.
Class \code{"\linkS4class{Intervals_virtual}"}, by class "Genome\_intervals", distance 3.
Class \code{"\linkS4class{matrix}"}, by class "Genome\_intervals", distance 4.
Class \code{"\linkS4class{array}"}, by class "Genome\_intervals", distance 5.
Class \code{"\linkS4class{structure}"}, by class "Genome\_intervals", distance 6.
Class \code{"\linkS4class{vector}"}, by class "Genome\_intervals", distance 7, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "Genome_intervals_stranded", to = "character")}: ... }
    \item{distance\_to\_nearest}{\code{signature(from = "Genome_intervals_stranded", to = "Genome_intervals_stranded")}: ... }
    \item{interval\_complement}{\code{signature(x = "Genome_intervals_stranded")}: ... }
    \item{interval\_intersection}{\code{signature(x = "Genome_intervals_stranded")}: ... }
    \item{interval\_overlap}{\code{signature(query = "Genome_intervals_stranded", target = "Genome_intervals_stranded")}: ... }
    \item{interval\_union}{\code{signature(x = "Genome_intervals_stranded")}: ... }
    \item{strand}{\code{signature(x = "Genome_intervals_stranded")}: ... }
    \item{strand<-}{\code{signature(x = "Genome_intervals_stranded")}: ... }
	 }
}


\seealso{
	\code{\linkS4class{Genome_intervals}} the parent class without strand.
}
\examples{

# The "Genome_intervals_stranded" class
j <- new(
	 "Genome_intervals_stranded",
	 matrix(
	 	c(1,2,  
	 	  3,5,
	 	  4,6,
	 	  8,9
	 	  ),
	 	byrow = TRUE,
                ncol = 2
	 ),
	 closed = matrix(
	 			c(
	 				FALSE, FALSE,
	 				TRUE, FALSE,
	 				TRUE, TRUE,
	 				TRUE, FALSE
	 			 ),
	 			byrow = TRUE,
     			ncol = 2
     			),
     annotation = data.frame(
     				seq_name = factor( c("chr01","chr01", "chr02","chr02") ),
	 				strand = factor( c("+", "+", "+", "-") ),
	 				inter_base = c(FALSE,FALSE,FALSE,TRUE)
	 				)
	 )

## print
print(j)

## size of each interval as count of included bases 
size(j)

## close intervals left and right (canonical representation)
close_intervals(j)
}

\keyword{classes}
