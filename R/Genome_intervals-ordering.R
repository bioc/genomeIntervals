### sorting ranking and ordering

## xtfrm() is a generic auxiliary function that produces a numeric vector which
##  will sort in the same order as 'x'.
## xtfrm() orders Genome_intervals objects according to
## 1) seq_name
## 2) [start & !inter-base < [start & inter-base < (start & !inter-base < (start & inter-base
## 3) stop) & !inter-base < stop) & inter-base < stop] & !inter-base < stop] & inter-base
## 4) strand (for Genome_intervals_stranded object)
##
## We get automatically order() which calls xtfrm() and implements rank() and sort() based on order() (code was taken from IRanges)

setMethod("xtfrm", "Genome_intervals", function(x){
			a = x[,1] + 1*(!closed(x)[,1]) + 0.5*inter_base(x)
			b = x[,2] - 1*(!closed(x)[,2]) + 0.5*inter_base(x)
			o = order(seqnames(x), a, b)
			rv = integer(length(o))
			rv[o] = seq_len(length(o))
			rv
		}
)

setMethod("xtfrm", "Genome_intervals_stranded", function(x){
			a = x[,1] + 1*(!closed(x)[,1]) + 0.5*inter_base(x)
			b = x[,2] - 1*(!closed(x)[,2]) + 0.5*inter_base(x)
			o = order(seqnames(x), a, b, strand(x))
			rv = integer(length(o))
			rv[o] = seq_len(length(o))
			rv
		}
)

### sort() relies on order() and on a "[" method for 'x'.
setMethod("sort", "Genome_intervals",
		function(x, decreasing=FALSE, ...)
		{
			x[order(x, decreasing=decreasing)]
		}
)

setMethod("rank", "Genome_intervals",
		function(x, na.last=TRUE, ties.method=c("average", "first", "last", "random", "max", "min"),...)
		{
			if (!missing(ties.method) && !identical(ties.method, "first"))
				stop("only 'ties.method=\"first\"' is supported when ranking GenomeIntervals")
			xtfrm(x)
		}
)
