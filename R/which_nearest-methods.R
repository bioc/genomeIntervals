setMethod("which_nearest",
          signature( "Genome_intervals", "Genome_intervals" ),
          function(from, to, check_valid = TRUE)
          {
            if ( check_valid && !( validObject(to) && validObject(from) ) )
              stop( "The 'to' and/or 'from' Genome_intervals are invalid." )
            
            ## result value, initiated as 3 empty objects (to be combined into data.frame later)
            n  <- nrow(from)
            rdist <- numeric(n)
            rnear <- vector("list", length=n)
            rover <- vector("list", length=n)
            
            ## adapt interval representation (basically moves to R with appropriate changes for inter_base)
            qints <- genomeIntervals:::intervalsForOverlap( to )
            tints <- genomeIntervals:::intervalsForOverlap( from)
            
            ## all unique seq_names
            seqlev <- unique( c( levels(seq_name(to)), levels(seq_name(from)) ) )
            ## loop over seq_names and call next method 
            for(s in seqlev){
                qi <- which(seq_name(to) == s)
                ti <- which(seq_name(from) == s)
                if (!any(ti) || !any(qi)) next
                wi <-  which_nearest(tints[ti], qints[qi],
                                     check_valid=check_valid)
                ## fill in respective result rows
                rdist[ti] <- wi$'distance_to_nearest'
                rnear[ti] <-
                  lapply(wi$'which_nearest', function(v) qi[v])
                rover[ti] <-
                  lapply(wi$'which_overlap', function(v) qi[v])
              }
            rv <- data.frame('distance_to_nearest'=rdist,
                             'which_nearest'=I(rnear),
                             'which_overlap'=I(rover),
                             row.names=rownames(from) )
            return(rv)
          }
) # setMethod("which_nearest") for two Genome_intervals objects


### TO DO if wanted: methods for stranded intervals
