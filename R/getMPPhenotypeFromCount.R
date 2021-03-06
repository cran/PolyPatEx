
##  Function getMPPhenotypeFromCount
##
## Given phenotype counts \code{nM}, \code{nP}, and \code{nMP},
## generate the corresponding maternal and progeny phenotypes.
##
## \code{getMPPhenotypeFromCount} generates the maternal and progeny
## phenotypes corresponding to the specified number of maternal
## alleles (\code{nM}), number of progeny alleles (\code{nP}), and
## number of alleles shared by progeny and mother (\code{nMP}).
##
## Maternal alleles (including those that also appear in the progeny)
## are represented as \code{M1}, \code{M2}, etc.  Progeny alleles
## that do not appear in its mother are represented as \code{P1},
## \code{P2}, etc.
##
## @title Generate maternal and progeny phenotypes
## @param nM integer: the number of alleles in the maternal phenotype
## @param nP integer: the number of alleles in the progeny's
## phenotype
## @param nMP integer: the number of alleles that appear in both
## maternal and progeny phenotypes
## @return A list with two elements:
##
## \describe{
##
## \item{\code{MPhenot}}{The maternal phenotype, as a character
## vector of alleles.}
##
## \item{\code{PPhenot}}{The progeny phenotype, as a character vector
## of alleles.}
##
## }
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## getMPPhenotypeFromCount(nM = 4, nP = 3, nMP = 1)
##
## }
##
getMPPhenotypeFromCount <- function(nM,nP,nMP) {
  ##
  if (nMP > min(nM,nP)) {  ## nMP impossibly large
    stop("\n nMP is larger than min(nM, nP)!\n\n")
  }
  if (nM == 0 || nP == 0 || nMP == 0) {
    ## Zero cases should be dealt with elsewhere...
    stop("\n One of nM, nP or nMP is zero!\n\n")
    ## Since this should only be called using parameters generated by
    ## getPossMPPhenotCounts(), this check is technically redundant.
  }
  MPhenot <- paste("M",1:nM,sep="")
  ## Sort order : M, then P
  PPhenot <- c(MPhenot[1:nMP],
               if (nMP < nP) paste("P",1:(nP-nMP),sep=""))
  return(list(MPhenot=MPhenot,PPhenot=PPhenot))
}
