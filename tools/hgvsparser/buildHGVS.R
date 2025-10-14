# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of hgvsParseR.
#
# hgvsParseR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hgvsParseR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hgvsParseR.  If not, see <https://www.gnu.org/licenses/>.

#' Genomic HGVS Builder 
#'
#' A constructor for a genomic-level HGVS builder object. The object contains a collection of functions
#' for building genomic HGVS strings.
#' 
#' The resulting object encapsulates the following functions:
#' \itemize{
#' 	\item{substitution(pos,ancestral,variant)} Genomic substitution variants. 
#'    pos = position (integer); ancestral = ancestral nucleotide [ACGT]; 
#'    variant = variant nucleotide [ACGT]
#' 	\item{deletion(start,stop)} Genomic deletion. start = start position (integer);
#'    stop = stop position (integer)
#' 	\item{inversion(start,stop)} Genomic inversion. start = start position (integer);
#'    stop = stop position (integer)
#' 	\item{duplication(start,stop)} Genomic duplication. start = start position (integer);
#'    stop = stop position (integer)
#' 	\item{insertion(start,variant)} Genomic insertion. start = position immediately preceeding 
#'    the insertion (integer); seq = inserted nucleotide sequence [ACGT]+
#' 	\item{delins(start,stop,variant)} Genomic deletion and insertion. start = start position (integer); 
#'    stop = stop position relative to the reference (integer); seq = inserted nucleotide sequence [ACGT]+
#'  \item{cis(...)} Multi-variant phased in cis. Parameters are genomic HGVS strings for the 
#'    corresponding single mutants
#'  \item{trans(...)} Multi-variant phased in trans. Parameters are genomic HGVS strings for the 
#'    corresponding single mutants
#'  \item{nophase(...)} Multi-variant with unknown phasing. Parameters are genomic HGVS strings for the 
#'    corresponding single mutants
#' }
#' 
#' @return A \code{hgvs.builder.g} object with functions for building genomic HGVS strings. 
#'   The individual functions return single-element character vectors containing these strings.
#' @keywords HGVS builder
#' @export
#' @examples
#' builder <- new.hgvs.builder.g()
#' string1 <- builder$substitution(123,"A","G")
#' string2 <- builder$delins(123,129,"ATTG")
#' string3 <- with(builder,cis(substitution(123,"A","C"),substitution(231,"G","A")))

new.hgvs.builder.g <- function() {

	substitution <- function(pos,ancestral,variant) {
		if (!is.numeric(pos) || pos < 1) stop("position must be a positive integer")
		if (!is.character(ancestral) || !(ancestral %in% c("A","C","G","T"))) stop("ancestral must be single nucleotide")
		if (!is.character(variant) || !(variant %in% c("A","C","G","T"))) stop("variant must be single nucleotide")
		paste0("g.",pos,ancestral,">",variant)
	}

	deletion <- function(start,stop) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (start > stop) stop("start must be upstream of stop")
		paste0("g.",start,"_",stop,"del")
	}

	inversion <- function(start,stop) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (start > stop) stop("start must be upstream of stop")
		paste0("g.",start,"_",stop,"inv")
	}

	duplication <- function(start,stop) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (start > stop) stop("start must be upstream of stop")
		paste0("g.",start,"_",stop,"dup")
	}

	insertion <- function(start,seq) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) stop("variant must be nucleotide sequence")
		paste0("g.",start,"_",start+1,"ins",seq)
	}

	delins <- function(start,stop,seq) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (start > stop) stop("start must be upstream of stop")
		if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) stop("variant must be nucleotide sequence")
		paste0("g.",start,"_",stop,"delins",seq)
	}

	cis <- function(...) {
		strings <- list(...)
		if (!all(sapply(strings,is.character))) stop("all arguments must be HGVS strings")
		strings <- unlist(strings)
		if (!all(substr(strings,1,2)=="g.")) stop("all arguments must be genomic HGVS strings")
		bodies <- substr(strings,3,nchar(strings))
		paste0("g.[",paste(bodies,collapse=";"),"]")
	}

	trans <- function(...) {
		strings <- list(...)
		if (!all(sapply(strings,is.character))) stop("all arguments must be HGVS strings")
		strings <- unlist(strings)
		if (!all(substr(strings,1,2)=="g.")) stop("all arguments must be genomic HGVS strings")
		bodies <- substr(strings,3,nchar(strings))
		paste0("g.[",paste(bodies,collapse="];["),"]")
	}

	nophase <- function(...) {
		strings <- list(...)
		if (!all(sapply(strings,is.character))) stop("all arguments must be HGVS strings")
		strings <- unlist(strings)
		if (!all(substr(strings,1,2)=="g.")) stop("all arguments must be genomic HGVS strings")
		bodies <- substr(strings,3,nchar(strings))
		paste0("g.[",paste(bodies,collapse="(;)"),"]")
	}

	return(structure(list(
		substitution=substitution,
		deletion=deletion,
		inversion=inversion,
		duplication=duplication,
		insertion=insertion,
		delins=delins,
		cis=cis,
		trans=trans,
		nophase=nophase
	),class="hgvs.builder.g"))
}

print.hgvs.builder.g <- function() {
	cat("Genomic HGVS string builder. Use $ operator to access functions.")
}



#' Coding Sequence HGVS Builder 
#'
#' A constructor for a CDS (=coding sequence) HGVS builder object. The object contains a collection of functions
#' for building CDS HGVS strings.
#' The resulting object encapsulates the following functions:
#' \itemize{
#' 	\item{substitution(pos,ancestral,variant,posOffset=0)} CDS substitution variants. 
#'    pos = position (integer); ancestral = ancestral nucleotide [ACGT]; 
#'    variant = variant nucleotide [ACGT]; posOffset = offset from the position when
#'    crossing exon-intron borders (integer, defaults to 0)
#' 	\item{deletion(start,stop,startOffset=0,stopOffset=0)} CDS deletion. start = start position (integer);
#'    stop = stop position (integer); startOffset = offset from the start position when
#'    crossing exon-intron borders (integer, defaults to 0); stopOffset = offset from the 
#'    stop position when crossing exon-intron borders (integer, defaults to 0)
#' 	\item{inversion(start,stop,startOffset=0,stopOffset=0)} CDS inversion. start = start position (integer);
#'    stop = stop position (integer); startOffset = offset from the start position when
#'    crossing exon-intron borders (integer, defaults to 0); stopOffset = offset from the 
#'    stop position when crossing exon-intron borders (integer, defaults to 0)
#' 	\item{duplication(start,stop,startOffset=0,stopOffset=0)} CDS duplication. start = start position (integer);
#'    stop = stop position (integer); startOffset = offset from the start position when
#'    crossing exon-intron borders (integer, defaults to 0); stopOffset = offset from the 
#'    stop position when crossing exon-intron borders (integer, defaults to 0)
#' 	\item{insertion(start,variant,startOffset=0)} CDS insertion. start = position immediately preceeding 
#'    the insertion (integer); seq = inserted nucleotide sequence [ACGT]+ ; 
#'    startOffset = offset from the start position when crossing exon-intron borders 
#'    (integer, defaults to 0)
#' 	\item{delins(start,stop,variant,startOffset=0,stopOffset=0)} CDS deletion and insertion. start = start position (integer); 
#'    stop = stop position relative to the reference (integer); 
#'    seq = inserted nucleotide sequence [ACGT]+ ; startOffset = offset from the start position when
#'    crossing exon-intron borders (integer, defaults to 0); stopOffset = offset from the 
#'    stop position when crossing exon-intron borders (integer, defaults to 0)
#'  \item{cis(...)} Multi-variant phased in cis. Parameters are coding HGVS strings for the 
#'    corresponding single mutants
#'  \item{trans(...)} Multi-variant phased in trans. Parameters are coding HGVS strings for the 
#'    corresponding single mutants
#'  \item{nophase(...)} Multi-variant with unknown phasing. Parameters are coding HGVS strings for the 
#'    corresponding single mutants
#' }
#' 
#' @return A \code{hgvs.builder.c} object with functions for building coding HGVS strings. 
#'   The individual functions return single-element character vectors containing these strings.
#' @keywords HGVS builder
#' @export
#' @examples
#' builder <- new.hgvs.builder.c()
#' string1 <- builder$substitution(123,"A","G",posOffset=2)
#' string2 <- builder$delins(123,129,"ATTG")
#' string3 <- with(builder,cis(substitution(123,"A","C"),substitution(231,"G","A")))

new.hgvs.builder.c <- function() {

	offsetStr <- function(offset) {
		if (offset==0) {
			""
		} else if (offset > 0) {
			paste0("+",offset)
		} else if (offset < 0) {
			as.character(offset)
		}
	}

	substitution <- function(pos,ancestral,variant,posOffset=0) {
		if (!is.numeric(pos) || pos < 1) stop("position must be a positive integer")
		if (!is.numeric(posOffset)) stop("offset must be an integer")
		if (!is.character(ancestral) || !(ancestral %in% c("A","C","G","T"))) stop("ancestral must be single nucleotide")
		if (!is.character(variant) || !(variant %in% c("A","C","G","T"))) stop("variant must be single nucleotide")
		paste0("c.",pos,offsetStr(posOffset),ancestral,">",variant)
	}

	deletion <- function(start,stop,startOffset=0,stopOffset=0) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (!is.numeric(startOffset)) stop("offset must be an integer")
		if (!is.numeric(stopOffset)) stop("offset must be an integer")
		if (start+startOffset > stop+stopOffset) stop("start must be before stop")
		paste0("c.",start,offsetStr(startOffset),"_",stop,offsetStr(stopOffset),"del")
	}

	inversion <- function(start,stop,startOffset=0,stopOffset=0) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (!is.numeric(startOffset)) stop("offset must be an integer")
		if (!is.numeric(stopOffset)) stop("offset must be an integer")
		if (start+startOffset >= stop+stopOffset) stop("start must be before stop")
		paste0("c.",start,offsetStr(startOffset),"_",stop,offsetStr(stopOffset),"inv")
	}

	duplication <- function(start,stop,startOffset=0,stopOffset=0) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (!is.numeric(startOffset)) stop("offset must be an integer")
		if (!is.numeric(stopOffset)) stop("offset must be an integer")
		if (start+startOffset > stop+stopOffset) stop("start must be before stop")
		paste0("c.",start,offsetStr(startOffset),"_",stop,offsetStr(stopOffset),"dup")
	}

	insertion <- function(start,seq,startOffset=0) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) {
			stop("variant must be nucleotide sequence")
		}
		if (!is.numeric(startOffset)) stop("offset must be an integer")
		stop <- if (startOffset != 0) start else start+1
		stopOffset <- if (startOffset != 0) startOffset+1 else startOffset
		paste0("c.",start,offsetStr(startOffset),"_",stop,offsetStr(stopOffset),"ins",seq)
	}

	delins <- function(start,stop,seq,startOffset=0,stopOffset=0) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (!is.numeric(startOffset)) stop("offset must be an integer")
		if (!is.numeric(stopOffset)) stop("offset must be an integer")
		if (start+startOffset > stop+stopOffset) stop("start must be before stop")
		if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) stop("variant must be nucleotide sequence")
		paste0("c.",start,offsetStr(startOffset),"_",stop,offsetStr(stopOffset),"delins",seq)
	}

	cis <- function(...) {
		strings <- list(...)
		if (!all(sapply(strings,is.character))) stop("all arguments must be HGVS strings")
		strings <- unlist(strings)
		if (!all(substr(strings,1,2)=="c.")) stop("all arguments must be coding HGVS strings")
		bodies <- substr(strings,3,nchar(strings))
		paste0("c.[",paste(bodies,collapse=";"),"]")
	}

	trans <- function(...) {
		strings <- list(...)
		if (!all(sapply(strings,is.character))) stop("all arguments must be HGVS strings")
		strings <- unlist(strings)
		if (!all(substr(strings,1,2)=="c.")) stop("all arguments must be coding HGVS strings")
		bodies <- substr(strings,3,nchar(strings))
		paste0("c.[",paste(bodies,collapse="];["),"]")
	}

	nophase <- function(...) {
		strings <- list(...)
		if (!all(sapply(strings,is.character))) stop("all arguments must be HGVS strings")
		strings <- unlist(strings)
		if (!all(substr(strings,1,2)=="c.")) stop("all arguments must be coding HGVS strings")
		bodies <- substr(strings,3,nchar(strings))
		paste0("c.[",paste(bodies,collapse="(;)"),"]")
	}

	return(structure(list(
		substitution=substitution,
		deletion=deletion,
		inversion=inversion,
		duplication=duplication,
		insertion=insertion,
		delins=delins,
		cis=cis,
		trans=trans,
		nophase=nophase
	),class="hgvs.builder.c"))
}

print.hgvs.builder.c <- function() {
	cat("Coding-sequence HGVS string builder. Use $ operator to access functions.")
}



#' Protein HGVS Builder 
#'
#' A constructor for a protein-level HGVS builder object. The object contains a collection of functions
#' for building protein HGVS strings.
#' 
#' The resulting object encapsulates the following functions:
#' \itemize{
#'  \item{synonymous()} A synonymous variant. No parameters required.
#'  \item{synonymous(pos,ancestral)} Unofficial (yet frequently used) version of synonymous variant syntax. 
#'    pos = position (integer); ancestral = ancestral amino acid in one-letter or three-letter code.
#' 	\item{substitution(pos,ancestral,variant)} AA substitution variants. 
#'    pos = position (integer); ancestral = ancestral amino acid in one-letter or three-letter code; 
#'    variant = variant amino acid in one-letter or three-letter code
#' 	\item{deletion(startPos,startAA,endPos,endAA)} AA deletion. startPos = start position (integer);
#'    startAA = start amino acid in one-letter or three-letter code;
#'    endPos = stop position (integer); endAA = start amino acid in one-letter or three-letter code
#' 	\item{duplication(startPos,startAA,endPos,endAA)} AA duplication. startPos = start position (integer);
#'    startAA = start amino acid in one-letter or three-letter code;
#'    endPos = stop position (integer); endAA = start amino acid in one-letter or three-letter code
#' 	\item{insertion(leftPos,leftAA,rightAA,seq)} AA insertion. leftPos = position immediately preceeding 
#'    the insertion (integer); leftAA = corresponding amino acid in one-letter or three-letter code;
#'    rightAA = amino acid to the right of the insertion, in one-letter or three-letter code;
#'    seq = inserted amino acid sequence, given as a character vector containing the individual
#'    one-letter or three-letter amino acid codes.
#' 	\item{delins(startPos,startAA,endPos,endAA,seq)} AA deletion and insertion. 
#'    startPos = start position (integer);
#'    startAA = start amino acid in one-letter or three-letter code;
#'    endPos = stop position (integer); endAA = start amino acid in one-letter or three-letter code;
#'    seq = inserted amino acid sequence, given as a character vector containing the individual
#'    one-letter or three-letter amino acid codes.
#'  \item{frameshift(startPos,startAA,variantAA=NA,newStop=NA)} Frameshift variant. 
#'    startPos = start position (integer);
#'    startAA = start amino acid in one-letter or three-letter code;
#'    variantAA = amino acid replacing the start position in the frameshift sequence, 
#'    given in one-letter or three-letter code, or \code{NA} to omit (default);
#'    newStop = the position of the nearest coding resulting from the frameshift, 
#'    or \code{NA} to omit (default).
#'  \item{cis(...)} Multi-variant phased in cis. Parameters are coding HGVS strings for the 
#'    corresponding single mutants. As phasing in trans would be nonsensical in a protein context,
#'    the \code{trans()} and \code{nophase()} methods are not provided here.
#' }
#' 
#' @return A \code{hgvs.builder.g} object with functions for building genomic HGVS strings. 
#'   The individual functions return single-element character vectors containing these strings.
#' @keywords HGVS builder
#' @export
#' @examples
#' builder <- new.hgvs.builder.g()
#' string1 <- builder$substitution(123,"R","K")
#' string2 <- builder$delins(123,"Arg",152,"Leu",c("Lys","Trp","Ser"))
#' string3 <- with(builder,cis(substitution(123,"R","K"),deletion(125,"S",152,"L")))


new.hgvs.builder.p <- function(aacode=c(1,3)) {

	aacode <- aacode[[1]]
	if (!is.numeric(aacode) && !(aacode %in% c(1,3))) {
		stop("Invalid aacode parameter, only 1 or 3 allowed!")
	}
	
	one2three <- c(A="Ala",C="Cys",D="Asp",E="Glu",F="Phe",G="Gly",H="His",
		I="Ile",K="Lys",L="Leu",M="Met",N="Asn",P="Pro",Q="Gln",R="Arg",
		S="Ser",T="Thr",V="Val",W="Trp",Y="Tyr",`*`="Ter")
	three2one <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",
		Gly="G",His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",
		Ser="S",Thr="T",Trp="W",Tyr="Y",Val="V",Ter="*")

	enforceCode <- function(aa) {
		if (aa %in% one2three) {
			if (aacode == 1) {
				three2one[[aa]]
			} else {
				aa
			}
		} else if (aa %in% three2one) {
			if (aacode == 1) {
				aa
			} else {
				one2three[[aa]]
			}
		} else {
			stop("Invalid AA code")
		}
	}

	synonymous <- function(pos=NULL,ancestral=NULL) {
		if (is.null(pos) || is.null(ancestral)) {
			return("p.=")
		}
		if (!is.numeric(pos) || pos < 1) stop("position must be a positive integer")
		if (!is.character(ancestral) || !(ancestral %in% c(one2three,three2one))) stop("ancestral must be single amimo acid")
		ancestral <- enforceCode(ancestral)
		paste0("p.",ancestral,pos,"=")
	}

	substitution <- function(pos,ancestral,variant) {
		if (!is.numeric(pos) || pos < 1) stop("position must be a positive integer")
		if (!is.character(ancestral) || !(ancestral %in% c(one2three,three2one))) stop("ancestral must be single amimo acid")
		if (!is.character(variant) || !(variant %in% c(one2three,three2one))) stop("variant must be single amino acid")
		ancestral <- enforceCode(ancestral)
		variant <- enforceCode(variant)
		paste0("p.",ancestral,pos,variant)
	}

	deletion <- function(startPos,startAA,endPos,endAA) {
		if (!is.numeric(startPos) || startPos < 1) stop("position must be a positive integer")
		if (!is.numeric(endPos) || endPos < 1) stop("position must be a positive integer")
		if (startPos > endPos) stop("start must be upstream of stop")
		if (!is.character(startAA) || !(startAA %in% c(one2three,three2one))) stop("startAA must be single amimo acid")
		if (!is.character(endAA) || !(endAA %in% c(one2three,three2one))) stop("endAA must be single amimo acid")
		startAA <- enforceCode(startAA)
		endAA <- enforceCode(endAA)
		if (startPos==endPos) {
			paste0("p.",startAA,startPos,"del")
		} else {
			paste0("p.",startAA,startPos,"_",endAA,endPos,"del")
		}
	}

	duplication <- function(startPos,startAA,endPos,endAA) {
		if (!is.numeric(startPos) || startPos < 1) stop("position must be a positive integer")
		if (!is.numeric(endPos) || endPos < 1) stop("position must be a positive integer")
		if (startPos >= endPos) stop("start must be upstream of stop")
		if (!is.character(startAA) || !(startAA %in% c(one2three,three2one))) 
			stop("startAA must be single amimo acid")
		if (!is.character(endAA) || !(endAA %in% c(one2three,three2one))) 
			stop("endAA must be single amimo acid")
		startAA <- enforceCode(startAA)
		endAA <- enforceCode(endAA)
		paste0("p.",startAA,startPos,"_",endAA,endPos,"dup")
	}

	insertion <- function(leftPos,leftAA,rightAA,seq) {
		if (!is.numeric(leftPos) || leftPos < 1) stop("position must be a positive integer")
		if (!is.character(leftAA) || !(leftAA %in% c(one2three,three2one))) 
			stop("leftAA must be single amimo acid")
		if (!is.character(rightAA) || !(rightAA %in% c(one2three,three2one))) 
			stop("rightAA must be single amimo acid")
		if (!is.character(seq) || !all(sapply(seq,function(x) x %in% c(one2three,three2one))))
			stop("seq must be a vector of amino acids")
		rightPos <- leftPos+1
		leftAA <- enforceCode(leftAA)
		rightAA <- enforceCode(rightAA)
		seq <- paste(sapply(seq,enforceCode),collapse="")
		paste0("p.",leftAA,leftPos,"_",rightAA,rightPos,"ins",seq)
	}

	delins <- function(startPos,startAA,endPos,endAA,seq) {
		if (!is.numeric(startPos) || startPos < 1) stop("position must be a positive integer")
		if (!is.numeric(endPos) || endPos < 1) stop("position must be a positive integer")
		if (startPos > endPos) stop("start must be upstream of stop")
		if (!is.character(startAA) || !(startAA %in% c(one2three,three2one))) 
			stop("startAA must be single amimo acid")
		if (!is.character(endAA) || !(endAA %in% c(one2three,three2one))) 
			stop("endAA must be single amimo acid")
		if (!is.character(seq) || !all(sapply(seq,function(x) x %in% c(one2three,three2one))))
			stop("seq must be a vector of amino acids")
		startAA <- enforceCode(startAA)
		endAA <- enforceCode(endAA)
		seq <- paste(sapply(seq,enforceCode),collapse="")
		paste0("p.",startAA,startPos,"_",endAA,endPos,"delins",seq)
	}

	frameshift <- function(startPos,startAA,variantAA=NA,newStop=NA) {
		if (!is.numeric(startPos) || startPos < 1) stop("position must be a positive integer")
		if (!is.na(newStop) && (!is.numeric(newStop) || newStop < 1)) stop("position must be a positive integer")
		if (!is.character(startAA) || !(startAA %in% c(one2three,three2one))) 
			stop("startAA must be single amimo acid or NA")
		if (!is.na(variantAA) && (!is.character(startAA) || !(startAA %in% c(one2three,three2one)))) 
			stop("variantAA must be single amimo acid or NA")
		startAA <- enforceCode(startAA)
		if (is.na(variantAA)) {
			variantAA <- ""
		} else {
			variantAA <- enforceCode(variantAA)
		}
		if (is.na(newStop)) {
			newStop <- "" 
		} else {
			newStop <- paste0("*",newStop)
		}
		paste0("p.",startAA,startPos,variantAA,"fs",newStop)
	}
	
	cis <- function(...) {
		strings <- list(...)
		if (!all(sapply(strings,is.character))) stop("all arguments must be HGVS strings")
		strings <- unlist(strings)
		if (!all(substr(strings,1,2)=="p.")) stop("all arguments must be protein HGVS strings")
		bodies <- substr(strings,3,nchar(strings))
		paste0("p.[",paste(bodies,collapse=";"),"]")
	}

	return(structure(list(
		synonymous=synonymous,
		substitution=substitution,
		deletion=deletion,
		duplication=duplication,
		insertion=insertion,
		delins=delins,
		frameshift=frameshift,
		cis=cis
	),class="hgvs.builder.p"))
}