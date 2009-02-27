# the size of a genome interval is in terms of its bases, i.e. computed on integers 'Z'
# inter_base intervals have size 0

setMethod(
        "size",
        signature( "Genome_intervals" ),
        function( x ) {
            rv <- callNextMethod(x)
            rv[ inter_base(x) ] <- 0 
            rv
        }
)



