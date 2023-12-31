\name{Genome_intervals-class}
\docType{class}
\alias{Genome_intervals-class}
\alias{[,Genome_intervals-method}
\alias{[,Genome_intervals,ANY-method}
\alias{[<-,Genome_intervals,ANY,missing,Genome_intervals-method}
\alias{[[,Genome_intervals-method}
\alias{[[,Genome_intervals,ANY-method}
\alias{[[,Genome_intervals,ANY,ANY-method}
\alias{[[<-,Genome_intervals-method}
\alias{[[<-,Genome_intervals,ANY-method}
\alias{[[<-,Genome_intervals,ANY,ANY-method}
\alias{$,Genome_intervals-method}
\alias{$<-,Genome_intervals-method}


\alias{annotation}
\alias{annotation<-}
\alias{annotation,Genome_intervals-method}
\alias{annotation<-,Genome_intervals-method}

\alias{inter_base}
\alias{inter_base<-}

\alias{coerce,Genome_intervals,Intervals_full-method}
\alias{coerce,Genome_intervals,character-method}
\alias{coerce,Genome_intervals,data.frame-method}

\alias{inter_base,Genome_intervals-method}
\alias{inter_base<-,Genome_intervals-method}

\alias{seqnames,Genome_intervals-method}
\alias{seqnames<-,Genome_intervals-method}
\alias{size,Genome_intervals-method}
\alias{show,Genome_intervals-method}
\alias{strand<-,Genome_intervals_stranded,ANY-method}
\alias{type<-,Genome_intervals-method}
\alias{which_nearest,Genome_intervals,Genome_intervals-method}

\alias{width,Genome_intervals-method}

\title{Class "Genome\_intervals" }

\description{ A set of genomic intervals without specified strand.
Genomic intervals are intervals over the integers with two further annotations: \code{seqnames} (a chromosome or more generally a sequence of origin)
and \code{inter_base} (\code{logical}) that states whether the interval is to be understood as an interval over bases
(such as coding-sequence) or inter-bases (such as restriction sites or insertion positions). }

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{See \code{\linkS4class{Intervals_full}} }
    \item{\code{annotation}:}{A \code{"data.frame"} with the same number of rows as \code{.Data}.
    	It has a column named \code{seq_name} that is a factor and does not contain missing values.
    	\code{seq_name} is used to represent the chromosome or more generally the sequence of origin of the intervals.
    	\code{annotation} has a column named \code{inter_base} that is logical and does not contain missing values.
    	\code{inter_base} is \code{FALSE} if the interval is to be understood as an interval over bases
(such as coding-sequence) and \code{TRUE} if it is over inter-bases (such as restriction site or an insertion position).
Like base intervals, inter-base interval are encoded over the integers. An inter-base at position \code{n} indicates the space between base \code{n} and \code{n+1}.
    	}
    \item{\code{closed}:}{See \code{\linkS4class{Intervals_full}} }
    \item{\code{type}:}{See \code{\linkS4class{Intervals_full} } }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Intervals_full}"}, directly.
Class \code{"\linkS4class{Intervals_virtual}"}, by class "Intervals\_full", distance 2.
Class \code{"\linkS4class{matrix}"}, by class "Intervals\_full", distance 3.
Class \code{"\linkS4class{array}"}, by class "Intervals\_full", distance 4.
Class \code{"\linkS4class{structure}"}, by class "Intervals\_full", distance 5.
Class \code{"\linkS4class{vector}"}, by class "Intervals\_full", distance 6, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "Genome_intervals")}: ... }
    \item{[[}{\code{signature(x = "Genome_intervals")}: ... }
    \item{[[<-}{\code{signature(x = "Genome_intervals")}: ... }
    \item{\$}{\code{signature(x = "Genome_intervals")}: ... }
    \item{\$<-}{\code{signature(x = "Genome_intervals")}: ... }

    \item{annotation}{\code{signature(object = "Genome_intervals")}: ... }
    \item{annotation<-}{\code{signature(object = "Genome_intervals")}: ... }

    \item{coerce}{\code{signature(from = "Genome_intervals", to = "Intervals_full")}: ... }
    \item{coerce}{\code{signature(from = "Genome_intervals", to = "character")}: ... }
    \item{coerce}{\code{signature(from = "Genome_intervals", to = "data.frame")}: ... }
    \item{distance\_to\_nearest}{\code{signature(from = "Genome_intervals", to = "Genome_intervals")}: ... }
    \item{inter\_base}{\code{signature(x = "Genome_intervals")}: ... }
    \item{inter\_base<-}{\code{signature(x = "Genome_intervals")}: ... }
    \item{interval\_complement}{\code{signature(x = "Genome_intervals")}: ... }
    \item{interval\_intersection}{\code{signature(x = "Genome_intervals")}: ... }
    \item{interval\_overlap}{\code{signature(from = "Genome_intervals", to = "Genome_intervals")}: ... }
    \item{interval\_union}{\code{signature(x = "Genome_intervals")}: ... }
    \item{seqnames}{\code{signature(x = "Genome_intervals")}: ... }
    \item{seqnames<-}{\code{signature(x = "Genome_intervals")}: ... }
    \item{size}{\code{signature(x = "Genome_intervals")}: ... }
    \item{type<-}{\code{signature(x = "Genome_intervals")}: ... }
    \item{which\_nearest}{
      For each interval in Set1, finds nearest (least distant)
      interval in Set2. Intervals on different chromsomes are
      never considered 'near' to each other.
      The returned value is a \code{data.frame} with the number
      of rows equal to the number of intervals in Set1. Each row
      specifies the distance to the nearest interval in Set2 (a
      0 means that the interval overlaps one or more intervals in Set2),
      and the indices of near and overlapping intervals in Set2.
      See \code{\linkS4class{Intervals_full}} for further details.}
    \item{width}{
      Returns the interval length as the number of bp covered (base interval) or
      spanned(inter-base interval). Similar to the \code{IRanges} package
      \code{\link[IRanges:width]{width}}} function.
    }
}

\note{
A \code{Genome_intervals} is a \code{"\linkS4class{Intervals_full}"} of type Z (i.e. a set of intervals over the integers).
The \code{annotation} slot can carry further columns that can serve as annotations.
}

\seealso{
	\code{\linkS4class{Genome_intervals_stranded}} for a derived class that allows stranded genomic intervals.
}
\examples{
# The "Genome_intervals" class

i <- new(
	 "Genome_intervals",
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
	 				TRUE, FALSE,
	 				TRUE, FALSE,
	 				TRUE, TRUE,
	 				TRUE, FALSE
	 			 ),
	 			byrow = TRUE,
     			ncol = 2
     			),
	 annotation = data.frame(
	 				seq_name = factor(c("chr01","chr01", "chr02","chr02")),
	 				inter_base = c(FALSE, FALSE, TRUE, TRUE)
	 				)
	 )

colnames(i) <- c( "start", "end" )

# print
print(i)

# size (number of bases per interval)
size(i)

## convert to a data.frame
as(i,"data.frame")

## simpler way to construct a Genome_intervals object:
G <- GenomeIntervals(start=c(1,3,4,5,10,8), end=c(5,5,6,8,11,9),
                     chromosome=rep(c("chr2","chrX","chr1"), each=2),
                     leftOpen=rep(c(FALSE, FALSE, TRUE), 2))
show(G)
}
\keyword{classes}
