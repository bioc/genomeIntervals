# complement of intervals by seq_name, strand and inter_base status
# just a wrapper for interval_complement

setMethod(
        "interval_complement",
        signature( "Genome_intervals" ),
        function( x ) {
            ##split by sequence and inter_base and call next method
            s = split(	x, list( inter_base(x), seqnames(x) ), drop=FALSE )

            fac.comb = expand.grid( levels(as.factor(inter_base(x))), levels(seqnames(x)) )

            interv.list = lapply( s , function(y) interval_complement( as(y,"Intervals_full") ) )
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
                                    levels = levels( seqnames(x) )
                            )
                    )
            )
        }
)


setMethod (
        "interval_complement",
        signature( "Genome_intervals_stranded" ),
        function( x ){
            ## split by strand and call next method
            s = split(x, strand(x), drop=TRUE )
            gi.list = lapply( s , function(y) interval_complement( as(y,"Genome_intervals") ) )
            gi = Reduce(c, gi.list)
            new(
                    "Genome_intervals_stranded",
                    gi@.Data,
                    closed = closed(gi),
                    annotation = data.frame(
                            seq_name = seqnames(gi),
                            inter_base = inter_base(gi),
                            strand= factor(
                                    rep(names(s), times = sapply( gi.list,nrow)),
                                    levels = levels( strand(x) )
                            )
                    )
            )
        }
)
