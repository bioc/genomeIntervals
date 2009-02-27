# benchmarking tests
# 
# on modifications, regenerate .Rout.save with:
# R --vanilla <./genomeIntervals/tests/benchmarking-tests.R >& ./genomeIntervals/tests/benchmarking-tests.Rout.save
# and use svn' diff to check for differences
#
# Author: julien.gagneur
###############################################################################

library(genomeIntervals)

options(warn = -1)

#---------------
# make data
#---------------
# toy examples
data("gen_ints")

# non-stranded versions
i0 = as(i, "Genome_intervals")
j0 = as(j, "Genome_intervals")
k0 = as(k, "Genome_intervals")

# empty intervals
e = k[1:3,]
e[,2] = e[,1]
closed(e) = FALSE

#---------------
# checks
#---------------
# distance to nearest
cat("distance_to_nearest\n")
print( distance_to_nearest(i,j) )		# x
print( distance_to_nearest(j,i) )		# x
print( distance_to_nearest(i,k) )		# x
print( distance_to_nearest(i,e) )		# x
print( distance_to_nearest(e,i) )		# x

print( distance_to_nearest(i0,j0) )		# x
print( distance_to_nearest(j0,i0) )		# x
cat("\n")

# interval overlap
cat("interval_overlap\n")
print( interval_overlap(i,k) )			# x
print( interval_overlap(i,e) )			# x
print( interval_overlap(i0,j0) )		# x
cat("\n")


# set operations
cat("interval_union\n")
print( close_intervals( interval_union(i) ) )		# x
print( close_intervals( interval_union(i,k) ) )		# x
print( close_intervals( interval_union(i0,k0) ) )	# x
print( close_intervals( interval_union(e) ) )		# x
cat("\n")

cat("interval_intersection\n")
print( close_intervals( interval_intersection(i,j) ) ) # x
print( close_intervals( interval_intersection(i,k) ) ) # x
print( close_intervals( interval_intersection(i,e) ) ) # x
cat("\n")

cat("interval_complement\n")
print( close_intervals( interval_complement(j) ) )	# x
print( close_intervals( interval_complement(j0) ) )	# x
print( close_intervals( interval_complement(k) ) ) # x
print( close_intervals( interval_complement(e) ) ) # x
cat("\n")

