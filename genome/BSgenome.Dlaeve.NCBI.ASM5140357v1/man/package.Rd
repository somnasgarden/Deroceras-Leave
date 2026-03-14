\name{BSgenome.Dlaeve.NCBI.ASM5140357v1}
\docType{package}

\alias{BSgenome.Dlaeve.NCBI.ASM5140357v1-package}
\alias{BSgenome.Dlaeve.NCBI.ASM5140357v1}
\alias{Dlaeve}

\title{Full genomic sequences for Deroceras laeve (NCBI assembly ASM5140357v1)}

\description{
  Full genomic sequences for Deroceras laeve as provided by NCBI (assembly ASM5140357v1, assembly accession GCA_051403575.1). The sequences are stored in DNAString objects.
}

\author{rlopezt <rafae@users.noreply.github.com>}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome}
          software package.

    \item The \code{\link[GenomeInfoDb]{seqinfo}} getter and
          \link[GenomeInfoDb]{Seqinfo} objects in the \pkg{GenomeInfoDb}
          package.

    \item The \code{\link[GenomeInfoDb]{seqlevelsStyle}} getter and
          setter in the \pkg{GenomeInfoDb} package.

    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.

    \item The \link[BSgenomeForge]{forgeBSgenomeDataPkgFromNCBI}
          function in the \pkg{BSgenomeForge} package that was used to
          create this package.
  }
}

\examples{
BSgenome.Dlaeve.NCBI.ASM5140357v1
bsg <- BSgenome.Dlaeve.NCBI.ASM5140357v1
head(seqlengths(bsg))
seqinfo(bsg)

## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
