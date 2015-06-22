# Unions of intervals by seq_name, strand and inter_base status
# just a wrapper for interval_union

setMethod(
        "interval_union",
        signature( "Genome_intervals" ),
        function( x, ... ) {

            args <- c( list(x), list(...) )

            ## drop non-core annotation
            core = lapply(args, core_annotated)

            ## combines all passed object into a single one
            x_comb <- do.call( c, core)

            ## then split by sequence and inter_base and call next method
            s = split(	x_comb, list( inter_base(x_comb), seqnames(x_comb) ), drop=FALSE )

            fac.comb = expand.grid( levels(as.factor(inter_base(x_comb))), levels(seqnames(x_comb)) )

            interv.list = lapply( s , function(y) interval_union( as(y,"Intervals_full") ) )
            nrows  = sapply(interv.list,nrow)

            interv = do.call(c, interv.list)

            new(
                    "Genome_intervals",
                    interv@.Data,
                    closed = closed(interv),
                    annotation = data.frame(
                            inter_base = rep( fac.comb[,1]=="TRUE", times = nrows),
                            seq_name = factor(
                                    rep( as.character(fac.comb[,2]), times = nrows),
                                    levels = levels( seqnames(x_comb) )
                            )
                    )
            )
        }
)


setMethod (
        "interval_union",
        signature( "Genome_intervals_stranded" ),
        function( x, ... ) {

            args <- c( list(x), list(...) )

            ## have args the same class?
            same_class <- all( sapply( args, is, "Genome_intervals_stranded" ) )
            if(!same_class)
                stop("All arguments should have the same class.")

            ## NAs in strand?
            strandNA <- sapply(args, function(y) any(is.na(strand(y))))
            if( any(strandNA) )
                stop("NA(s) present in the strand of at least one argument.")

            ## drop non-core annotation
            core = lapply(args, core_annotated)

            ## combines all passed object into a single one
            x_comb <- do.call( c, core)

            ## then split by strand and call next method
            s = split(x_comb, strand(x_comb), drop=TRUE )
            gi.list = lapply( s , function(y) interval_union( as(y,"Genome_intervals") ) )
            gi = do.call(c, gi.list)
            new(
                    "Genome_intervals_stranded",
                    gi@.Data,
                    closed = closed(gi),
                    annotation = data.frame(
                            seq_name = seqnames(gi),
                            inter_base = inter_base(gi),
                            strand= factor(
                                    rep(names(s), times = sapply( gi.list,nrow)),
                                    levels = levels( strand(x_comb) )
                            )
                    )
            )
        }
)
