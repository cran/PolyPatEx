
##  Function removeMismatches
##
## Remove mismatching cases between mothers and their progeny.
##
## \code{removeMismatches} calls \code{\link{getProgenyStatusGenot}}
## or \code{\link{getProgenyStatusPhenot}} (as appropriate) to
## identify cases where a progeny's mother is incapable of providing
## a gamete compatible with the progeny's allele set.
##
## Where only one locus in a progeny-mother pair is a mismatch
## (\sQuote{single locus mismatch}), the offending locus is set to
## have no alleles (a \sQuote{missing locus}) in the progeny.  This
## locus will therefore be ignored in subsequent analyses.
##
## Where more than one locus is a mismatch(\sQuote{multi-locus
## mismatch}), the progeny is removed entirely from the dataset.
##
## \code{removeMismatches} will report details of such mismatches to
## the R console, to help the user to track down and check their
## origins, and also allow them to see whether the problem may lie
## with the mother, rather than in her progeny.
##
## @title Remove mismatches between mother and progeny
## @param adata data frame: an allele dataset.
## @return An allele dataset as a data frame, with identified
##              mismatches removed (see details)
## @author Alexander Zwart (alec.zwart at csiro.au)
##
removeMismatches <- function(adata) {
  ##
  checkForValidPPEDataset(adata)
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dataType <- attr(adata,"dataType")
  dioecious <- attr(adata,"dioecious")
  selfCompatible <- attr(adata,"selfCompatible")
  ##
  if (dataType=="genotype") {
    progenyStatusTable <- getProgenyStatusGenot(adata)$progenyStatusTable
  } else {  ## Phenotypic data
    progenyStatusTable <- getProgenyStatusPhenot(adata)$progenyStatusTable
  }
  ##
  locusMismatches <- progenyStatusTable == "MP.noMatch"
  singleLocusMismatches <- apply(locusMismatches,1,sum) == 1
  multiLocusMismatches <- apply(locusMismatches,1,sum) > 1
  ##
  if (any(multiLocusMismatches)) {
    multiMismatchProgeny <- names(multiLocusMismatches)[multiLocusMismatches]
    cat("\n Multi-locus (i.e., 2 or more) mismatches between mother and")
    cat("\n progeny have been found.  The relevant progeny will be")
    cat("\n removed from the dataset.  The progeny and their mismatching")
    cat("\n loci are: \n\n")
    ML.Ptab <-  cbind(ProgenyId = multiMismatchProgeny,
                      Loci = apply(locusMismatches[multiLocusMismatches,,drop=FALSE],
                        1,
                        function(vv){
                          paste(which(vv),collapse=", ")
                        })
                      )
    rownames(ML.Ptab) <- NULL
    print(ML.Ptab,quote=FALSE)
    ##Tally the progeny from each mother that are to be removed.
    t1 <- with(adata,table(mother[id %in% multiMismatchProgeny]))
    ##The total number of progeny for each each mother
    t2 <- table(stripNAs(adata$mother))
    ML.Mtab <- data.frame(Mother=names(t1),
                          num.ML.mismatches=as.vector(t1),
                          total.progeny = as.vector(t2[match(names(t1),
                            names(t2))]))
    cat("\n The following table provides a summary for the mothers of the")
    cat("\n progeny above.  The table shows, for each affected mother,")
    cat("\n the number of progeny in which such multilocus mismatches")
    cat("\n occur, and the total number of progeny for that mother. This")
    cat("\n may help identify possible problems with a mother's allele set.\n\n")
    print(ML.Mtab,quote=FALSE)
    ##
    ##Remove progeny with multi-locus mismatches from the dataset.
    adata <- adata[!(adata$id %in% multiMismatchProgeny),]
    ##
  }
  ##
  if (any(singleLocusMismatches)) {
    ##
    singleMismatchProgeny <- names(singleLocusMismatches)[singleLocusMismatches]
    cat("\n\n Single-locus mismatches between mother and progeny have")
    cat("\n been found.  The relevant progeny-locus combinations will")
    cat("\n be set to contain no alleles ('allele set missing'). The")
    cat("\n progeny and their mismatching loci are: \n\n")
    cc <-  cbind(ProgenyId = singleMismatchProgeny,
                 Locus = apply(locusMismatches[singleLocusMismatches,,drop=FALSE],
                   1,
                   function(vv){
                     paste(which(vv),collapse=", ")
                     })
                 )
    rownames(cc) <- NULL
    cc <- as.data.frame(cc)
    cc$Locus <- as.numeric(as.character(cc$Locus))
    print(cc,quote=FALSE)
    cat("\n")
    ##
    ##Tally the progeny from each mother that are to be 'repaired'.
    t1 <- with(adata,table(mother[id %in% singleMismatchProgeny]))
    ##The total number of progeny for each each mother
    t2 <- table(stripNAs(adata$mother))
    SL.Mtab <- data.frame(Mother=names(t1),
                          num.SL.mismatches=as.vector(t1),
                            total.progeny = as.vector(t2[match(names(t1),
                              names(t2))]))
    cat("\n The following table summarises the mothers having progeny with")
    cat("\n single locus mismatches, the number of progeny per mother in")
    cat("\n which such mismatches occur, and the total number of progeny")
    cat("\n per mother.  This may help identify possible problems with a")
    cat("\n mother's allele set.\n\n")
    print(SL.Mtab,quote=FALSE)
    cat("\n")
    ##
    ##Set any single-locus mismatches to be all missing.
    for (thisProgeny in as.character(cc$ProgenyId)) {
      locusRange <- (3 + dioecious) + 1:ploidy +
        (with(cc,Locus[ProgenyId==thisProgeny])-1)*ploidy
      is.na(adata[thisProgeny, locusRange]) <- TRUE
    }
    cat("\n Note that single-locus mismatches will subsequently be")
    cat("\n recorded as 'P.missing', since such loci have had all")
    cat("\n progeny alleles removed from the dataset\n\n")
  }
  ##
  return(adata)
}
