## adapt intervals for genome interval overlap operation

intervalsForOverlap = function(x){
    rv = as(x, "Intervals_full")
    
    ## close the non-empty intervals and set the empty ones to NA
    e = empty(rv)
    if( any(e) ) warning("Empty intervals encountered, set to NA.", call.=FALSE)
    
    rv[!e,] = close_intervals(rv[!e,])
    
    # set NA to empty intervals
    rv[e,] = NA
    
    ## set intervals type to R
    type(rv) <- "R"
    
    ## represent inter_base intervals as right open left open intervals (requires adding one to the right bound))
    rv[inter_base(x)][,2] <- rv[inter_base(x)][,2] + 1
    closed(rv)[inter_base(x),] <- FALSE
    
    return(rv)
}


setMethod(
        "interval_overlap",
        signature( "Genome_intervals", "Genome_intervals" ),
        function(
                target,
                query,
                check_valid = TRUE
        ){
            if ( check_valid && !( validObject(query) && validObject(target) ) )
                stop( "The 'query' and/or 'target' Genome_intervals are invalid." )
            
            ## result value, initiated as empty list
            rv = vector( mode="list", length=nrow(target) )
            
            ## adapt interval representation (basically moves to R with appropriate changes for inter_base)
            qints <- intervalsForOverlap( query )
            tints <- intervalsForOverlap( target)
            
            ## all unique seq_names
            seqlev = unique( c( levels(seq_name(query)), levels(seq_name(target)) ) )
            
            ## loop over seq_names and call next method 
            for(s in seqlev){
                qi = which( seq_name(query) == s)
                ti = seq_name(target) == s
                
                rv[ti] = lapply(
                        interval_overlap(tints[ti], qints[qi], check_valid=check_valid),
                        function(v) qi[v]
                )
            }
            
            return(rv)
        }
)


setMethod(
        "interval_overlap",
        signature( "Genome_intervals_stranded", "Genome_intervals_stranded" ),
        function(
                target,
                query,
                check_valid = TRUE
        ){
            if ( check_valid && !( validObject(query) && validObject(target) ) )
                stop( "The 'query' and/or 'target' Genome_intervals are invalid." )
            
            if( any( is.na(strand(query)) ) || any( is.na(strand(target)) ) )
                stop("NA(s) present in the strand of 'query' or 'target'.")
            ## result value, initiated as empty list
            rv = vector( mode="list", length=nrow(target) )
            strdlev = levels( strand(target) )
            ## loop over strands and call next method
            nextMethod = getMethod("interval_overlap", signature( "Genome_intervals",  "Genome_intervals" ) )
            for(s in strdlev){
                qi = which(strand(query) == s)
                ti = strand(target) == s
                
                rv[ti] =  lapply(
                        nextMethod( target[ti], query[qi] ),
                        function(v) qi[v]
                )
            }
            return(rv)
        }
)
