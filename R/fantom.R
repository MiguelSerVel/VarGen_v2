# ---- FANTOM ----

#' @title Generate data.frame from FANTOM5 promoters file
#' @description Prepare FANTOM5 promoters file for
#'  \code{\link{get_fantom5_variants}}
#'
#' @param promoters_file the "hg38_liftover+new_CAGE_peaks_phase1and2_ann.txt" 
#' file from FANTOM5.
#' The file can be downloaded here:
#' https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_liftover+new_CAGE_peaks_phase1and2_ann.txt.gz
#' @return a data.frame containing the file information with the following columns:
#' \itemize{
#'   \item genome (either hg19 or hg38)
#'   \item chrom (chromosome)
#'   \item chromStart (start position of the promoter site)
#'   \item chromEnd (end position of the enhancer site)
#'   \item name (promoter ID)
#'   \item symbol (hgnc symbol)
#'   }
prepare_fantom_promoters <- function(promoters_file) {
  # Read promoters table
  promoters <- utils::read.delim(promoters_file, header = T,
                                 comment.char = "#", stringsAsFactors = FALSE)
  
  descriptions <- promoters$short_description
  
  pattern <- "(.*)::(.*):(\\d+)\\.\\.(\\d+),([+-]);(.*)"
  promoters <- strcapture(pattern,
                          promoters[,1],
                          data.frame(genome=character(),
                                     chrom=character(),
                                     chromStart=integer(),
                                     chromEnd=integer(),
                                     strand=character(),
                                     name=character()))
  
  # Get HGNC symbols
  promoters$symbol <- sub(".*@", "", descriptions)

  return(promoters)
}


#' @title Generate data.frame from FANTOM5 enhancers file
#' @description Prepare FANTOM5 promoters file for
#'  \code{\link{get_fantom5_variants}}
#'
#' @param promoters_file the "F5.hg38.enhancers.bed" file from FANTOM5.
#' The file can be downloaded here:
#' https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz
#' @param score_threshold the minimum number of CAGE peaks detected for that enhancer.
#' This allows to filter out low atictivy enhancers (default: 5)
#' @return a data.frame containing the file information with the following columns:
#' \itemize{
#'   \item chrom (chromosome)
#'   \item chromStart (start position of the enhancer site)
#'   \item chromEnd (end position of the enhancer site)
#'   \item name (string with chromosome, start and end positions)
#'   \item score ( = correlation*1000)
#'   \item strand (+, -, or . (for unknown or not strand-specific enhancers))
#'   \item thichStart (start of the "thick" region used in UCSC Genome Browser)
#'   \item thickEnd (end of the "thick" region used in UCSC Genome Browser)
#'   \item itemRgb (RGB colour code)
#'   \item blockCount (number of blocks)
#'   \item blockSizes (comma-separated list of block sizes)
#'   \item chromStarts (comma-separated list of start offsets)
#'   }
prepare_fantom_enhancers <- function(enhancers_file, score_threshold = 5){
  # Get enhancers from file
  enhancers <- utils::read.delim(enhancers_file, header = F,
                                 stringsAsFactors = FALSE)
  
  colnames(enhancers) <- c("chrom", "chromStart", "chromEnd", "name", "score",
                           "strand", "thickStart", "thickEnd", "itemRgb",
                           "blockCount", "blockSizes", "chromStarts")
  
  # Get the enhancers with a decent level of correlations & significance
  # Correlation is Pearson as a z-score:
  # "(pearson corr - (mean of random motifs)) / std(pearson of random motifs)"
  # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0560-6
  # A z-score greater than 0 represents an element greater than the mean, this
  # means "more correlation than random motifs"
  # score = correlation*1000
  enhancers <- enhancers[(enhancers$score >= score_threshold),]
  
  return(enhancers)
}


#' @title Get variants on the enhancers of the list of genes given as input
#' @description FANTOM5 is used to get the enhancers of the genes, then the variants
#' located on the enhancers are fetched with \code{\link{get_variants_from_locations}}
#' It is also possible to specify a correlation threshold to limit the number
#' of association to output (default: 0.25). The enhancer tss association file
#' is based on hg19, so liftOver is needed.
#'
#' @param fantom_df the output of \code{\link{prepare_fantom}}
#' @param omim_genes output from \code{\link{get_omim_genes}}
#' @param enhancer_gap the length in bases to look for enhancers from a 
#' given chromosomal position for promoters (default: 100000)
#' @param hg19ToHg38.over.chain the chain file to liftOver locations from
#' hg19 to hg38.
#' @param verbose if true, will print progress information (default: FALSE)
#'
#' @return a data.frame with information about the variants located on the enhancers.
#' The data.frame will contain the following columns:
#' \itemize{
#'   \item chr (chromosome)
#'   \item pos (position of the variant)
#'   \item rsid (variant ID)
#'   \item ensembl_gene_id ("gene id" of the gene associated with the variant)
#'   \item hgnc_symbol ("hgnc symbol" of the gene associated with the variant)
#'   \item source (here the value will be "fantom5")
#' }
#
# gene_mart <- connect_to_gene_ensembl()
# DM1_genes <- get_omim_genes(omim_ids = "222100", gene_mart = gene_mart)
#
# vargen_install(install_dir = "./vargen_data/")
# fantom_df <- prepare_fantom("./vargen_data/enhancer_tss_associations.bed")
# get_fantom5_variants(fantom_df, DM1_genes, 0.25, "hg19ToHg38.over.chain")
get_fantom5_variants <- function(promoters_df, enhancers_df, omim_genes,
                                 hg19ToHg38.over.chain, enhancer_gap = 100000,
                                 verbose = FALSE) {
  fantom_variants <- data.frame()
  list.variants <- vector('list', nrow(omim_genes))

  for(gene in 1:nrow(omim_genes)){
    
    # Get only promoters from the genes of interest
    promoters <- promoters_df[(promoters_df$symbol == omim_genes[gene, "hgnc_symbol"]),]
    
    if(nrow(promoters) != 0) {
      # Get promoters with hg19 build and hg38 build
      promoters_19 <- promoters[(promoters$genome == "hg19"),]
      promoters_38 <- promoters[(promoters$genome == "hg38"),]
      
      # Convert each data.frame to granges object
      promoters_19_gr <- GenomicRanges::GRanges(seqnames = promoters_19$chrom,
                                                ranges = IRanges::IRanges(promoters_19$chromStart,
                                                                          promoters_19$chromEnd))
      promoters_38_gr <- GenomicRanges::GRanges(seqnames = promoters_38$chrom,
                                                ranges = IRanges::IRanges(promoters_38$chromStart,
                                                                          promoters_38$chromEnd))
      
      # Perform liftover ono hg19 build and join both objects
      promoters_19_gr <- unlist(rtracklayer::liftOver(promoters_19_gr, 
                                                      rtracklayer::import.chain(hg19ToHg38.over.chain)))
      promoters_gr <- c(promoters_19_gr, promoters_38_gr)
      
      # Obtain genomic ranges object for enhancers
      enhancers_gr <- GenomicRanges::GRanges(seqnames = enhancers_df$chrom,
                                             ranges = IRanges::IRanges(enhancers_df$chromStart,
                                                                       enhancers_df$chromEnd))
      
      # Take the position of enhancers that are 100kb from the promoters
      nearby_pairs <- GenomicRanges::findOverlaps(promoters_gr, enhancers_gr, 
                                                  maxgap = enhancer_gap)
      
      enhancers_gr <- enhancers_gr[S4Vectors::subjectHits(nearby_pairs)]
      
      # Obtain locations for promoters and enhancers
      promoters_locs <- paste0(GenomeInfoDb::seqnames(promoters_gr), ":",
                               BiocGenerics::start(promoters_gr)-1, ":",
                               BiocGenerics::end(promoters_gr))
      promoters_locs <- sub("^chr", "", promoters_locs)
      
      enhancers_locs <- paste0(GenomeInfoDb::seqnames(enhancers_gr), ":",
                               BiocGenerics::start(enhancers_gr)-1, ":",
                               BiocGenerics::end(enhancers_gr))
      enhancers_locs <- sub("^chr", "", enhancers_locs)
      
      fantom_locs <- c(promoters_locs, enhancers_locs)
      
      # Filter fantom_locs
      fantom_locs <- fantom_locs[fantom_locs != "::"]
      
      # Get variants from all locations found
      variants_locs <- get_variants_from_locations(fantom_locs,
                                                   verbose = TRUE)
      variants_locs <- unique(variants_locs)
      
    
      if(length(variants_locs) != 0){
        enhancer_variants <- cbind(ensembl_gene_id = omim_genes[gene, "ensembl_gene_id"],
                                   hgnc_symbol = omim_genes[gene, "hgnc_symbol"],
                                   variants_locs)

        enhancer_variants_df <- format_output(chr = unlist(enhancer_variants$seq_region_name),
                                              pos =  unlist(enhancer_variants$start),
                                              rsid = unlist(enhancer_variants$id),
                                              ensembl_gene_id = unique(enhancer_variants$ensembl_gene_id),
                                              hgnc_symbol = unique(enhancer_variants$hgnc_symbol))

        list.variants[[gene]] <- enhancer_variants_df
      }
    }
  }

  # Removing the NULL elements of the list if exists
  list.variants <- list.variants[!sapply(list.variants, is.null)]
  # Then concatenate the list into a data.frame
  fantom_variants <- do.call('rbind', list.variants)
  if(length(fantom_variants) != 0){
    fantom_variants$source <- "fantom5"
  }

  return(fantom_variants)
}
