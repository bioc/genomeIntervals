##' Coercion methods of the genomeIntervals package
##'
##' \describe{
##' \item{coerce}{ This method allows to coerce a
##' \code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals
##' object}} into a \code{\link[IRanges:RangedData-class]{RangedData object}} or
##' \code{\linkS4class{GRangesList}} object.}}
##'
##' @aliases coerce,Genome_intervals,RangedData-method
##' coerce,Genome_intervals,GRangesList-method
##' coerce,Genome_intervals,GRanges-method
##' as,Genome_intervals,RangedData-method
##' as,Genome_intervals,GRangesList-method
##' as as,Genome_intervals-method
##' @name genomeIntervals coercion methods
##' @rdname Genome_intervals-coercion-methods
##' @param from An object of class \code{\linkS4class{Genome_intervals}}
##' @param to a character string; either RangedData or GRangesList
##' @usage \S4method{as}{Genome_intervals}(from,to)
##' @return \describe{
##' \item{coerce}{ A \code{\linkS4class{RangedData}} or
##' \code{\linkS4class{GRangesList}}
##' containing the result of the coercion.  }}
##' @examples
##' \dontrun{
##' annot<-readGff3(system.file("extdata",
##'                             "Dmel-mRNA-exon-r5.52.gff3",
##'                             package="RnaSeqTutorial")
##' gAnnot<-as(annot,"RangedData") type(annot)
##' }
##' @author Nicolas Delhomme
##' @seealso
##' \itemize{
##' \item \code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals
##' object}}
##' \item \code{\link[genomeIntervals:readGff3]{readGff3 function}}
##' }

## convert a genome intervals into a RangedData
## TODO find a way to check and keep unexpected slots.
## or at least warn they would be ignored
setAs(from="Genome_intervals",to="RangedData",def=function(from){
  universe="intervals"

  # first check
  if (!is(from, "Genome_intervals")){stop("'from' must be a Genome_intervals object")}

  # create the ranges
  ranges<-IRanges(start=from[,1],end=from[,2])

  ## create the space
  ## drop the original levels
  space = as.character(seqnames(from))
  names(ranges) <- seq(along=ranges)

  # create the values
  slotToKeep<-na.omit(match(c("strand","exon","feature","intron","transcript","transcript.name","gene","gene.name"),names(from@annotation)))
  values <- as(data.frame(from@annotation[slotToKeep]),"DataFrame")
  elementMetadata(values)<-as(data.frame(labelDescription=names(from@annotation)[slotToKeep],row.names=names(from@annotation)[slotToKeep]),"DataFrame")

  # important before splitting
  rownames(values) <- names(ranges)

  # now split properly
  if (length(unique(space)) > 1) {
    ranges <- split(ranges, space)
    values <- split(values, space)
  } else {
    ranges <- IRangesList(ranges)
    values <- SplitDataFrameList(values, compress = TRUE)
    names(ranges) <- unique(space)
    names(values) <- names(ranges)
  }

  if (!is.null(universe) && !S4Vectors::isSingleString(universe))
    stop("'universe' must be a single string")
  universe(ranges) <- universe
  return(new(
      "RangedData",
      ranges = ranges,
      values = values))
})

## coerce into GRanges/List
setAs(from="Genome_intervals",to="GRanges",def=function(from){

  ## first check
  if (!is(from, "Genome_intervals")){stop("'from' must be a Genome_intervals object")}

  ## get all possible gff attributes
  ## to be quick
  ## we expect the types to have the same annotation
  ## so we fetch the first instance of each type and
  ## get the gff attribute names
  ## and use these to get the gffAttributes
  ## and convert that into a df
  mat <- do.call(cbind,lapply(unique(unlist(lapply(
    parseGffAttributes(from[match(unique(as.character(type(from))),type(from))]),
    names))),function(attr,from){
      getGffAttribute(from,attr)
    },from))

  ## create the object
  return(GRanges(ranges=IRanges(
    start=from[,1],
    end=from[,2]),
    seqnames=seqnames(from),
    strand=strand(from),
    cbind(
      DataFrame(apply(
        annotation(from)[!colnames(annotation(from))
                         %in% c("seq_name","strand","gffAttributes")],2,Rle)),
      DataFrame(apply(mat,2,Rle)))
  ))
})


## coerce into
setAs(from="Genome_intervals",to="GRangesList",def=function(from){
  ## create the object
  return(split(as(from,"GRanges"),seqnames(from)))
})
