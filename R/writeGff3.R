## write Gff3
setGeneric(name="writeGff3",
           def=function(object,file=character(1)){
             standardGeneric("writeGff3")
           })

setMethod(f="writeGff3",
          signature="Genome_intervals",
          definition=function(object,file=character(1)){
            ## convert and pass to the next function
            writeGff3(as(object,"data.frame"),file)
          })

setMethod(f="writeGff3",
          signature="data.frame",
          definition=function(object,file=character(1)){

            ## check the file
            stopifnot(file != "")
            stopifnot(file.exists(dirname(file)))

            ## check that we have 9 columns
            if(ncol(object)!=9){
              stop("Your object is not valid to be written as a gff3 as it does not have 9 columns.")
            }

            ## and the proper columns
            if(!all(colnames(object) %in% c("seqname","source","feature","start","end","score","strand","frame","attribute"))){
              stop(paste("Your object is not valid to be written as a gff3 as it does not the 9 expected columns.",
                         "These are: seqname, source, feature, start, end, score, strand, frame, attribute",sep="\n"))
            }

            ## TODO check the type in addition if we really want to enfore gff3...

            ## make sure they are in the right order
            object <- object[,match(c("seqname","source","feature","start","end","score","strand","frame","attribute"),colnames(object))]

            ## write
            write("##gff-version 3",file=file)
            write.table(format(object,scientific = FALSE),
                        append=TRUE,quote=FALSE,sep="\t",
                        row.names=FALSE,col.names=FALSE,
                        file=file)
          })

