setMethod(
		"show",
		signature( "Genome_intervals" ),
		function( object ) {
			cat(
					"Object of class ",
					class( object ),
					"\n",
					sum( !inter_base(object) ),
					" base interval",
					ifelse( sum( !inter_base(object) ) == 1, "", "s" ),
					" and ",
					sum( inter_base(object) ),
					" inter-base interval",
					ifelse( sum( inter_base(object) ) == 1, "", "s" ),
					"(*)",
					":\n",
					sep = ""
			)
			
			## start and end part of the object
			nshow = 4
			if(nrow(object)>2*nshow){
				nrowStart <- min (nshow , nrow(object) )
				startObject <- object[1:nrowStart]
				
				nrowEnd <- min (nshow , nrow(object) - nrowStart)
				endObject <- object[(nrow(object)-nrowEnd+1):nrow(object)] 
			}else{
				startObject <- object
				nrowStart <- nrow(object)
				nrowEnd <- 0
			}
			
			## show intervals
			ints <- as( startObject, "character")
			if ( !is.null( rownames( startObject ) ) ) {
				fmt <- sprintf( "%%%is", max( nchar( rownames( startObject ) ) ) )
				ints <- paste( sprintf( fmt, rownames( startObject ) ), ints )
			}
			cat( ints, sep = "\n" )
			
			if(nrowEnd > 0 ){
				cat("...", nrow(object) - nrowStart - nrowEnd, "other intervals...\n" )
				
				ints <- as( endObject, "character")
				if ( !is.null( rownames( endObject ) ) ) {
					fmt <- sprintf( "%%%is", max( nchar( rownames( endObject ) ) ) )
					ints <- paste( sprintf( fmt, rownames( endObject ) ), ints )
				}
				cat( ints, sep = "\n" )
			}

			## show annotations
			cat( "\nannotation:\n")
			show(annotation(startObject))
			if(nrowEnd > 0 ){
				cat("...", nrow(object) - nrowStart - nrowEnd, "other intervals...\n" )
				show(annotation(endObject))
			}
		}
)
