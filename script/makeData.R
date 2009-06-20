# makeData.R
# 
# generate data objects for package genomeIntervals
#
# Author: julien.gagneur
###############################################################################


library(genomeIntervals)

i <- new(
        "Genome_intervals_stranded",
        matrix(
                c(
                        1,2,  
                        2,5,
                        11,12,
                        8,9,
                        5,10,
                        4,12,
                        2,6
                ),
                byrow = TRUE,
                ncol = 2
        ),
        closed = matrix(
                c(
                        TRUE, TRUE,
                        FALSE, FALSE,
                        TRUE, FALSE,
                        TRUE, FALSE,
                        TRUE, TRUE,
                        TRUE, TRUE,
                        TRUE, FALSE
                ),
                byrow = TRUE,
                ncol = 2
        ),
        annotation = data.frame(
                seq_name = c("chr01","chr01","chr02","chr02", "chr02", "chr02", "chr03"),
                inter_base = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
                strand = c("+", "+", "-", "-", "-","+", "+")
        )
)


j <- new(
        "Genome_intervals_stranded",
        matrix(
                c(1,2,  
                        5,10,
                        4,6,
                        12,15,
                        8,9
                ),
                byrow = TRUE,
                ncol = 2
        ),
        closed = matrix(
                c(
                        TRUE, FALSE,
                        FALSE, TRUE,
                        TRUE, FALSE,
                        TRUE, TRUE,
                        TRUE, FALSE
                ),
                byrow = TRUE,
                ncol = 2
        ),
        annotation = data.frame(
                seq_name = c("chr01","chr01", "chr01","chr02","chr02"),
                inter_base = c(FALSE,FALSE,FALSE,FALSE,FALSE),
                strand = c("+", "+", "-", "+", "-")
        )
)

rownames(i) = paste("i.gene", 1:nrow(i), sep=".")


rownames(j) = paste("j.gene", 1:nrow(j), sep=".")

## example with inter_base
k = j
inter_base(k) = c(TRUE,FALSE,FALSE,FALSE,TRUE)
rownames(k)[!inter_base(k)] = paste("k.gene", 1:sum( !inter_base(k) ), sep =".")
rownames(k)[inter_base(k)] = paste("k.site", 1:sum(inter_base(k)), sep =".")

save(list=ls(), file = "data/gen_ints.rda"  )


new(
    "Genome_intervals_stranded",
    matrix(c(1, 2, 2, 5), ncol = 2),
    closed = matrix(c(TRUE, TRUE, TRUE, FALSE), ncol = 2),
    annotation = data.frame(
            seq_name = c("chr01","chr02"),
            inter_base = c(FALSE,FALSE),
            strand = factor( c("+", "+"), levels=c("+", "-") )
    )
)


