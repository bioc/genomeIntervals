# consistency test
#
# Author: julien.gagneur
###############################################################################


library( genomeIntervals )

options(warn = -1)

#---------
# settings
#---------
# size of random objects
n = 1e+3
# chrom length
l = 1e+3
# nb of chroms
k = 3

#---------
# data generation
#---------

randGenint = function(n,l,k){
    m = matrix( sample(l, 2*n, replace=TRUE), nc = 2 )
    m = cbind( apply( m, 1, min), apply( m, 1, max ) )

    cl = matrix( sample( c(FALSE,TRUE), 2*n, replace =TRUE), nc=2 )
    new(
        "Genome_intervals_stranded",
        m,
        closed = cl,
        annotation = data.frame(
                seq_name = paste("chr", sample(k, n, replace=TRUE) ),
                inter_base = sample( c(FALSE,TRUE), n, replace =TRUE),
                strand = sample(c("-", "+"), n, replace =TRUE)
        )
    )
}

i = randGenint(n,l,k)
j = randGenint(n,l,k)
i0 = as(i, "Genome_intervals")
j0 = as(j, "Genome_intervals")

#---------
# checks
#---------

# distances from i to j
dn = distance_to_nearest(i,j)

# distance is NA or >=0
if( any( !is.na(dn) & dn < 0) ) stop("negative distance.")

# distance == 0 if and only if the interval overlaps another one:
io = interval_overlap(i,j)
if( any( ( sapply(io, length) >0 )  != (!is.na(dn) & dn ==0) ) )
    stop("The property 'distance == 0 if and only if the interval overlaps another one' is not followed for at least one instance.")

# same test for not stranded objects
dn0 = distance_to_nearest(i0,j0)
if( any( !is.na(dn0) & dn0 < 0) ) stop("negative distance.")

io = interval_overlap(i0,j0)
if( any( ( sapply(io, length) >0 )  != (!is.na(dn0) & dn0 ==0) ) )
    stop("The property 'distance == 0 if and only if the interval overlaps another one' is not followed for at least one instance.")

# unstranded distance <= stranded distance
delta = dn - dn0
if( any(!is.na(delta)  & delta < 0) ) stop("some unstranded distance larger than a stranded one.")

# intersection with complement is empty
stopifnot( nrow( interval_intersection(i, interval_complement(i) ) ) == 0 )

# distance of union with complement is 1
# test must be done for not inter-base (or inter-base) independently
a = interval_union(i[!inter_base(i),] )
b = interval_complement(i[!inter_base(i),])

if(!(all.equal( distance_to_nearest( a, b ), rep(1, nrow(a) )  ) ) )
    stop("distance of union with complement is not 1.")

# width is reported consistently
# they should all be 4 in length (we alternate the open/closed state of the intervals
# pairwise)
gi <- GenomeIntervals(start=c(6,6,5,5),
                      end=c(10,9,10,9),
                      chromosome=rep("chr1",4),
                      leftOpen = c(FALSE,FALSE,TRUE,TRUE),
                      rightOpen=c(TRUE,FALSE,TRUE,FALSE))

stopifnot(all(width(gi)==4))
