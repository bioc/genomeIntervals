### To update the DESCRIPTION
# installed.packages()[c("GenomeInfoDb","GenomicRanges","IRanges","S4Vectors"),"Version"]
### To update the man and NAMESPACE
# roxygenize(".",roclets=c('rd', 'collate', 'namespace'),clean=TRUE)
###==========================
## To define the NAMESPACE
###==========================
#### IMPORT
#### Generic import
##' @import methods
#### Classes
##' @importClassesFrom intervals Intervals_virtual Intervals_virtual_or_numeric
##'  Intervals_full
#### S4 Methods
##' @importMethodsFrom BiocGenerics annotation "annotation<-" Reduce
##' strand "strand<-" width
##' @importMethodsFrom GenomeInfoDb seqnames "seqnames<-"
##' @importMethodsFrom intervals "[" close_intervals closed "closed<-"
##' empty head close_intervals open_intervals interval_union interval_complement
##' interval_intersection interval_overlap distance_to_nearest which_nearest
##' size tail type "type<-"
##' @importMethodsFrom IRanges "universe<-"
##' @importMethodsFrom S4Vectors "elementMetadata<-" Rle
#### S3 Methods
##' @importFrom GenomicRanges GRanges GRangesList
##' @importFrom IRanges IRanges IRangesList SplitDataFrameList
##' @importFrom S4Vectors DataFrame
#### EXPORT
#### Classes
##' @exportClass Genome_intervals Genome_intervals_stranded
#### S4 Methods
##' @exportMethod annotation annotation<- c seq_name "seq_name<-" seqnames
##' "seqnames<-" strand "strand<-" inter_base "inter_base<-" "type<-" "[" "[["
##' "$" coerce size show rank sort xtfrm core_annotated interval_union
##' interval_complement interval_intersection interval_overlap
##' distance_to_nearest which_nearest width readBasePairFeaturesGff3
##' readZeroLengthFeaturesGff3 readGff3 writeGff3
#### S3 Methods
##' @export parseGffAttributes getGffAttribute GenomeIntervals
NULL

###==========================
## To detail the deprecation
###==========================
##' The following function have been deprecated:
##' \itemize{
##' \item \code{seq_name}
##' \item \code{seq_name<-}
##' }
##'
##' \itemize{
##' \item The \code{seq_name} and \code{seq_name<-} accessor functions have been
##' replaced by the more generic \code{seqnames} and \code{seqnames<-} accessor
##' functions, respectively.
##' }
##' @aliases seq_name seq_name<-
##' seq_name,Genome_intervals-method
##' seq_name<-,Genome_intervals-method
##' @name Genome_intervals deprecated functions
##' @rdname Genome_intervals-deprecated
NULL
