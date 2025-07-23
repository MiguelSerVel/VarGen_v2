# ---- FANTOM ----

#' @title Generate data.frame from FANTOM5 promoters file
#' @description Prepare FANTOM5 promoters file for
#'  \code{\link{get_fantom5_enhancers_from_hgnc}}
#'
#' @param promoters_file the "hg38_liftover+new_CAGE_peaks_phase1and2_ann.txt" 
#' file from FANTOM5.
#' The file can be downloaded here:
#' https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_liftover+new_CAGE_peaks_phase1and2_ann.txt.gz
#' @return a data.frame containing the file information
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
#'  \code{\link{get_fantom5_enhancers_from_hgnc}}
#'
#' @param promoters_file the "F5.hg38.enhancers.bed" file from FANTOM5.
#' The file can be downloaded here:
#' https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz
#' @return a data.frame containing the file information
prepare_fantom_enhancers <- function(enhancers_file, corr_threshold = 0.25){
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
  score_threshold <- corr_threshold*1000
  enhancers <- enhancers[(enhancers$score >= score_threshold),]
  
  return(enhancers)
}


#' @title Get the enhancers associated to certain genes from FANTOM5
#' @description from the FANTOM5 dataset, get the enhancers associated to the
#' genes using the HGNC symbols, the association must pass the correlation threshold given
#' by the user. Is used internally by \code{\link{get_fantom5_variants}}
#'
#' @param fantom_df the output of \code{\link{prepare_fantom}}
#' @param hgnc_symbols vector of HUGO ids for the genes of interest
#' @param corr_threshold the minimum correlation (z-score) to consider a
#' enhancer/gene assocation valid (default: 0.25).
#' A z-score greater than 0 represents an element greater than the mean, this means
#' that this association has more correlation than random motifs.
#' @return a subset of the FANTOM5 data.frame containing the information about the
#' enhancers. The data.frame contains the following columns
#' \itemize{
#'   \item chr (chromosome)
#'   \item start (start of the enhancers)
#'   \item end (end of the enhancers)
#'   \item symbol (HGNC symbol of the gene associated to the enhancer)
#'   \item corr (correlation z-score from FANTOM5)
#'   \item fdr (False Discovery Rate)
#' }
#'

# vargen_install(install_dir = "./vargen_data/")
# fantom_df <- prepare_fantom("./vargen_data/enhancer_tss_associations.bed")
# get_fantom5_enhancers_from_hgnc(fantom_df = fantom_df,
#                                 hgnc_symbols = "FOXP3",
#                                 corr_threshold = 0.30)
get_fantom5_enhancers_from_hgnc <- function(fantom_df, hgnc_symbols,
                                            corr_threshold = 0.25){
  # Get the enhancers with a decent level of correlations & significance
  # Correlation is Pearson as a z-score:
  # "(pearson corr - (mean of random motifs)) / std(pearson of random motifs)"
  # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0560-6
  # A z-score greater than 0 represents an element greater than the mean, this
  # means "more correlation than random motifs"
  enhancers_corr <- unique(subset(fantom_df,
                                  fantom_df$corr >= corr_threshold,# & fdr < fdr_threshold,
                                  select = c("chr", "start", "end", "symbol", "corr", "fdr"))
  )
  # Get list of enhancer related to list of genes
  return(enhancers_corr[(enhancers_corr$symbol %in% hgnc_symbols),])
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
#' @param corr_threshold the minimum correlation (z-score) to consider a
#' enhancer/gene association valid (default: 0.25). A z-score greater than 0
#' represents an element greater than the mean, this means that this association
#'  has more correlation than random motifs.
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
get_fantom5_variants <- function(fantom_df, omim_genes, corr_threshold = 0.25,
                                 hg19ToHg38.over.chain, verbose = FALSE) {
  fantom_variants <- data.frame()
  list.variants <- vector('list', nrow(omim_genes))

  for(gene in 1:nrow(omim_genes)){
    enhancers_df <- get_fantom5_enhancers_from_hgnc(fantom_df = fantom_df,
                                                    hgnc_symbols = omim_genes[gene, "hgnc_symbol"],
                                                    corr_threshold = corr_threshold)
    if(nrow(enhancers_df) != 0) {
      enhancers_df <- GenomicRanges::makeGRangesFromDataFrame(enhancers_df,
                                                              keep.extra.columns = TRUE)
      enhancers_df <- unlist(rtracklayer::liftOver(enhancers_df, rtracklayer::import.chain(hg19ToHg38.over.chain)))
      fantom_locs <- paste0(GenomeInfoDb::seqnames(enhancers_df), ":",
                            BiocGenerics::start(enhancers_df)-1, ":",
                            BiocGenerics::end(enhancers_df))
      fantom_locs <- sub("^chr", "", fantom_locs)


      variants_loc <- get_variants_from_locations(fantom_locs,
                                                  verbose = verbose)
      if(length(variants_loc) != 0){
        enhancer_variants <- cbind(ensembl_gene_id = omim_genes[gene, "ensembl_gene_id"],
                                   hgnc_symbol = omim_genes[gene, "hgnc_symbol"],
                                   variants_loc)

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
