### ======
### Import
### ======
## Gff3
setGeneric(name = "readGff3",
           def = function(file,
                          isRightOpen = FALSE,
                          quiet = FALSE){
             standardGeneric("readGff3")
             })

## Base Pair features
setGeneric(name = "readBasePairFeaturesGff3",
           def = function(file,
                          quiet = FALSE){
             standardGeneric("readBasePairFeaturesGff3")
           })

## Zero length features
setGeneric(name = "readZeroLengthFeaturesGff3",
           def = function(file,
                          quiet = FALSE){
             standardGeneric("readZeroLengthFeaturesGff3")
           })

### ======
### Export
### ======
setGeneric(name="writeGff3",
           def=function(object,file=character(1)){
             standardGeneric("writeGff3")
           })
