# the core (minimal) annotation of Genome_intervals object
# useful before calls to c()
# 
# Author: julien.gagneur
###############################################################################

setGeneric( "core_annotated", function(x) standardGeneric( "core_annotated" ) )

setMethod(
        "core_annotated",
        signature( "Genome_intervals" ),
        function( x ) {
            rv = x
            annotation(rv) <- annotation(x)[, c("seq_name", "inter_base")] 
            rv
        }
)

setMethod(
        "core_annotated",
        signature( "Genome_intervals_stranded" ),
        function( x ) {
            rv = x
            annotation(rv) <- annotation(x)[, c("seq_name", "inter_base", "strand")] 
            rv
        }
)

