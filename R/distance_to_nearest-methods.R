## adapt intervals for distance computation

intervalsForDistance = function(x){

    rv = as(x, "Intervals_full")

    ## close the non-empty intervals and set the empty ones to NA
    e = empty(rv)
    if( any(e) ) warning("Empty intervals encountered, set to NA.", call.=FALSE)

    rv[!e,] = close_intervals(rv[!e,])

    # set NA to empty intervals
    rv[e,] = NA

    ## set intervals type to R
    type(rv) <- "R"

    ## shift inter_base intervals by 0.5 to the right (makes sense only after moving to R)
    rv[inter_base(x)] <- rv[inter_base(x)] + 0.5

    return(rv)
}



setMethod(
        "distance_to_nearest",
        signature( "Genome_intervals", "Genome_intervals" ),
        function( from, to ) {
            ## return value
            rv = vector( mode="numeric", length = nrow(from) )

            ## adapt interval representation (basically moves to R with appropriate changes for inter_base)
            ## empty intervals become NA
            fints <- intervalsForDistance( from )
            tints <- intervalsForDistance( to )

            seqlev = levels( seqnames(from) )

            ## loop over seq_names and call next method
            for(s in seqlev){
                fi = seqnames(from) == s

                rv[fi] = distance_to_nearest(
                        fints[fi],
                        tints[seqnames(to) == s]
                )
            }
            return(rv)
        }
)


## if both are stranded, then take strands in account
## user must cast to unstranded if this is not wished.

setMethod(
        "distance_to_nearest",
        signature( "Genome_intervals_stranded", "Genome_intervals_stranded" ),
        function( from, to ) {
            if( any( is.na(strand(to)) ) )
                stop("NA(s) present in the strand of 'to'.")

            ## return value
            rv = vector( mode="numeric", length = nrow(from) )
            strdlev = levels( strand(from) )
            ## loop over strands and call next method
            nextMethod = getMethod("distance_to_nearest", signature( "Genome_intervals",  "Genome_intervals" ) )
            for(s in strdlev){
                fi = strand(from) == s
                rv[fi] =  nextMethod( from[fi], to[strand(to) == s] )
            }
            return(rv)
        }
)
