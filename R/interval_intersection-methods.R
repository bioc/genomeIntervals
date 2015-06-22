# common combination of factors
# idenitifed as alle entris in the first data.frame of factor combination
# that have at least one entry matching in each of the factor combination dataframes

commonCombination = function(df.list){
    common = sapply(
            1:nrow(df.list[[1]]),
            function(i){
                all(
                        sapply(
                                df.list,
                                function(df){
                                    any(
                                            sapply(
                                                    1:nrow(df),
                                                    function(j)
                                                        all( df[j,] == df.list[[1]][i,] )
                                            )
                                    )
                                }
                        )
                )
            }
    )
    return(df.list[[1]][common,])
}

## interval intersection has to be redefined (although complement and union are already there)
## because the input objects may not have the same seq_name and strand levels.
## In such case, the complement of an object missing one seq_name/strand level may be incomplete
## with respect to another one.

setMethod(
        "interval_intersection",
        signature( "Genome_intervals" ),
        function( x, ... ) {
            ## work by sequences and inter_base, restricting to the set of common sequences
            args <- c( list(x), list(...) )

            ## get unified seq_name levels
            seqlev = sort( unique(unlist(lapply(args, function(x) levels(seqnames(x) ) ) ) ) )

            ## combination of factors per arg
            fac.comb.list <- lapply(
                    args,
                    function(y){
                        rv = unique(annotation(y)[, c("seq_name", "inter_base")] )
                        rv$seq_name = as.character(rv$seq_name)
                        rv
                    }
            )

            fac.comb <- commonCombination( fac.comb.list)

            ## do the intersection by seq_names and inter_base
            interv.list = lapply(
                    1:nrow(fac.comb),
                    function(i){
                        ## extract the intervals
                        ints = lapply(
                                args,
                                function(y){
                                    as( y[ (inter_base(y)==fac.comb$inter_base[i]) &
                                             (seqnames(y)==fac.comb$seq_name[i]) ], "Intervals_full" )}
                        )
                        ## apply the intervals intersection
                        do.call( interval_intersection, ints)
                    }
            )
            nrows = lapply(interv.list, nrow)

            ## combine intervals
            interv = do.call(c, interv.list)

            ## make returned object
            new(
                    "Genome_intervals",
                    interv@.Data,
                    closed = closed(interv),
                    annotation = data.frame(
                            inter_base = rep( fac.comb$inter_base, times = nrows),
                            seq_name = factor( rep( fac.comb$seq_name, times = nrows), levels=seqlev )
                    )
            )
        }
)


setMethod(
        "interval_intersection",
        signature( "Genome_intervals_stranded" ),
        function( x, ... ) {

            ## work by sequences and inter_base, restricting to the set of common sequences
            args <- c( list(x), list(...) )

            same_class <- all( sapply( args, is, "Genome_intervals_stranded" ) )
            if(!same_class)
                stop("All arguments should have the same class.")

            strandNA <- sapply(args, function(y) any(is.na(strand(y))))
            if( any(strandNA) )
                stop("NA(s) present in the strand of at least one argument.")

            ## get unified seq_name levels
            seqlev = sort( unique(unlist(lapply(args, function(x) levels(seqnames(x))))) )

            ## get unified strd lev
            strdlev = sort( unique(unlist(lapply(args, function(x) levels(strand(x))))) )
            if(length(strdlev)!=2)
                stop("The arguments do not have the same levels for strand.")

            ## combination of factors per arg
            fac.comb.list <- lapply(
                    args,
                    function(y){
                        rv = unique(annotation(y)[, c("seq_name", "strand", "inter_base")] )
                        rv$seq_name = as.character(rv$seq_name)
                        rv$strand = as.character(rv$strand)
                        rv
                    }
            )

            fac.comb <- commonCombination( fac.comb.list)

            ## do the intersection by seq_names and inter_base
            interv.list = lapply(
                    1:nrow(fac.comb),
                    function(i){
                        ## extract the intervals
                        ints = lapply(
                                args,
                                function(y)
                                    as( y[
                                      (inter_base(y)==fac.comb$inter_base[i]) &
                                        (seqnames(y)==fac.comb$seq_name[i]) &
                                        (strand(y)==fac.comb$strand[i])
                                      ], "Intervals_full" )
                        )
                        ## apply the intervals intersection
                        do.call( interval_intersection, ints)
                    }
            )
            nrows = lapply(interv.list, nrow)

            ## combine intervals
            interv = do.call(c, interv.list)

            ## make returned object
            new(
                    "Genome_intervals_stranded",
                    interv@.Data,
                    closed = closed(interv),
                    annotation = data.frame(
                            inter_base = rep( fac.comb$inter_base, times = nrows),
                            seq_name = factor( rep( fac.comb$seq_name, times = nrows), levels =seqlev ),
                            strand = factor(rep(fac.comb$strand, times = nrows), levels =strdlev)
                    )
            )
        }
)
