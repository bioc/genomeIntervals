# the fullShow methods is a show method for Genome_intervals objects.
# It is called by test scripts so that the full output of the tests can be checked.  
# 
# Author: gagneur
###############################################################################

fullShow <- function(object ) {
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
	ints <- as( object, "character")
	if ( !is.null( rownames( object ) ) ) {
		fmt <- sprintf( "%%%is", max( nchar( rownames( object ) ) ) )
		ints <- paste( sprintf( fmt, rownames( object ) ), ints )
	}
	cat( ints, sep = "\n" )
	cat( "annotation:\n")
	show( annotation(object) )
}



