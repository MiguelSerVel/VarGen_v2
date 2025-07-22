# ---- VCF ----

#' @title studies the content of a vcf and finds specific genetic variants
#' @description loads a VCF file and looks for variants from a data.frame of 
#' genetic variants ids that can be annotaed or not. If not, these variants are
#' then annotated. If a GWAS file is provided, GWAS information related to the 
#' found variants is added to the data.frame. The user can also plot the variants 
#' with the function vargen_visualisation, highlighting those variants that have
#' been found in the VCF file.
#'
#' @param vcf_file a string with the path to a VCF file 
#' @param rsid_df a data.frame with genetic variants ids
#' @param outdir a string with the directory to store the images generated
#' @param min_qual the minimum quality required to filter VCF variants
#' @param memory_limit the number of RAM GB set for DuckDB to operate (default: 2)
#' @param max_headlines the maximum number of lines expected to have the VCF file
#' to search for the names of the columns (default: 1000)
#' @param df_annotated a boolean to specify if the "rsid_df" is annotated to avoid
#' re-annotating the data.frame (default: FALSE)
#' @param plot_vars a boolean to specify if the user wants to plot the variants
#' found
#' @param verbose if true, will print progress information (default: FALSE)
#'
#' @return a data.frame with the annotated variants from VCF file that are in the
#' "rsid_df"
#'
#' @examples
#' vargen_install("./vargen_data/")
#' 
#' DM1_simple <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid_ids = "222100",
#'                               fantom_corr = 0.25, outdir = "./", verbose = TRUE)
#'    
#' vcf_vars <- study_vcf("./path_to_vcf_file.vcf.gz", DM1_simple, outdir = "./",
#'                       plot_vars = TRUE, verbose = TRUE)                           
#' @export
study_vcf <- function(vcf_file, rsid_df, gwas_file, outdir = "./", min_qual = 20,
                      memory_limit = 2, max_headlines = 1000, df_annotated = FALSE,
                      plot_vars = FALSE, verbose = FALSE){
  
  # Create outdir if it does not exist
  if (!file.exists(outdir)){
    if(verbose) print(paste0("Creating folder '", outdir, "'"))
    dir.create(outdir)
  }
  
  
  #-----------------------------------------------------------------------------
  # Data frame extraction and annotation
  #-----------------------------------------------------------------------------
  
  # Load vcf file dataframe 
  if(verbose) print("Loading data from VCF file...")
  vcf_df <- load_vcf_file(vcf_file, memory_limit = memory_limit,
                          max_headlines = max_headlines)
  
  # Obtain dataframe with shared rsids between vcf file and rsid dataframe
  if(verbose) print("Comparing vcf file data with rsid dataframe...")
  shared_df <- compare_vcf(vcf_df = vcf_df, rsid_df = rsid_df,
                           memory_limit = memory_limit, verbose = verbose)
  
  # Drop rows with less quality than the minimum
  shared_df <- shared_df[shared_df$QUAL >= min_qual, ]
  
  if(!df_annotated){
    # Annotate dataframe with shared rsids
    if(verbose) print("Annotating dataframe with shared rsids...")
    shared_df <- annotate_dataframe(shared_df, verbose = verbose)
  } 
  
  # Filter dataframe by reference and alterantive bases
  shared_df <- shared_df[(shared_df$REF == shared_df$vcf_ref &
                            shared_df$ALT == shared_df$vcf_alt),]
  
  # Delete unnecessary columns
  shared_df$vcf_ref <- NULL
  shared_df$vcf_alt <- NULL
  
  shared_df <- unique(shared_df)
  
  # Create a column with the first word for the snpeff impact
  shared_df$impact_most_severe <- sapply(strsplit(shared_df$snpeff_ann, ";\\s*"),
                                         function(x) x[1])
  
  
  #-----------------------------------------------------------------------------
  # GWAS file annotation
  #-----------------------------------------------------------------------------
  
  if(!missing(gwas_file)){
    if(file.exists(gwas_file)){
      
      # Load gwas dataframe
      if(verbose) print("Loading data from GWAS file...")
      gwas_df <- load_gwas_file(gwas_file, memory_limit = memory_limit)
      
      # Obtain dataframe with shared rsids between already found rsids and rsids 
      # GWAS file
      if(verbose) print("Comparing GWAS file data with selected rsids...")
      intersect_df <- compare_gwas(gwas_df = gwas_df, rsid_df = shared_df,
                                   memory_limit = memory_limit, verbose = verbose)
      
      # Append new information to shared_df. GWAS rows for rsids not found in 
      # GWAS file will be NA
      shared_df <- merge(shared_df, intersect_df, 
                         by = intersect(names(shared_df), names(intersect_df)), 
                         all.x = TRUE)
      
    } else{
      warning("GWAS file does not exist. Skipping GWAS annotation.")
    }
  }
  
  
  #-----------------------------------------------------------------------------
  # Plotting
  #-----------------------------------------------------------------------------
  
  if(plot_vars){
    # Use vargen visualisation to plot variants found in vcf compated to variants
    # in master list
    
    # Check if rsid_df is already annotated
    if(!df_annotated){
      if(verbose) print("Generating plots. Variant annotation may take a while.")
      # Annotate original list of RSIDs
      rsid_df_ann <- annotate_dataframe(rsid_df, verbose = verbose)
    } else{
      if(verbose) print("Generating plots.")
      rsid_df_ann <- rsid_df
    }
    
    # Highlight rsids found in vcf file
    marked_rsids <- shared_df$rsid
    marked_rsids <- marked_rsids[!is.na(marked_rsids)]
    
    # Create plots
    vargen_visualisation(rsid_df_ann, outdir = outdir, rsid_highlight = marked_rsids,
                         device = "png", verbose = verbose)
  }
  
  if(verbose) print("Generating summary plot...")
  
  # Plot proportion of variants found depending on source
  rsid_df_small <- rsid_df[, c("rsid", "source")]
  shared_df_ann_small <- shared_df[, c("rsid", "source")]
  
  rsid_df_small$dataframe <- "Master list"
  shared_df_ann_small$dataframe <- "VCF file"
  
  all_variants <- rbind(rsid_df_small, shared_df_ann_small)
  
  # Count number of variants per source and gene
  count_df <- as.data.frame(table(all_variants$source, all_variants$dataframe))
  colnames(count_df) <- c("source", "dataframe", "frequency")
  
  # Remove zero counts
  count_df <- count_df[count_df$frequency > 0, ]
  
  # Plot
  count_plot <- ggplot2::ggplot(count_df, 
                                ggplot2::aes(x = source, y = frequency, 
                                             fill = dataframe)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_text(ggplot2::aes(label = frequency), 
                       position = ggplot2::position_dodge(width = 0.9), 
                       vjust = -0.5, size = 3) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +  
    ggplot2::labs(title = "Number of variants per source",
                  x = "Source", y = "Frequency")
  
  name = paste0(outdir, "/variant_freq_summary.jpeg")
  jpeg(name, width = 2000, height = 1200, res = 200)
  print(count_plot)  
  dev.off()
  
  
  #-----------------------------------------------------------------------------
  # Return annotated dataframe with shared rsids
  #-----------------------------------------------------------------------------
  
  return(shared_df)
}


#' @title Loads a VCF file into the R environment
#' @description Loads a VCF file into a data.frame in the R environment using
#' fast loading with DuckDB
#'
#' @param vcf_file a string with the path to the VCF file
#' @param memory_limit the number of RAM GB set for DuckDB to operate (default: 2)
#' @param max_headlines the maximum number of lines expected to have the VCF file
#' to search for the names of the columns (default: 1000)
#' @return a data.frame with the information stored in the VCF file
#'
#' @examples
#' vcf_df <- load_vcf_file("./path_to_vcf_file,vcf.gz")
#' @export
load_vcf_file <- function(vcf_file, memory_limit = 2, max_headlines = 1000){
  
  if (!file.exists(vcf_file)){
    stop("Couldn't find the file.")
  }
  
  # Path to vcf file
  file_name = vcf_file
  
  # Path to DuckDB file
  db_path <- "./duckdb_database.duckdb"
  
  # Open duckDB connection
  con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = db_path)
  
  # Set RAM usage limit to specified number of GB
  memory_query <- sprintf("SET memory_limit='%sGB'", memory_limit)
  DBI::dbExecute(con, memory_query)
  
  # Use read_csv_auto with custom options
  query <- sprintf("
    SELECT * FROM read_csv_auto('%s', delim='\t', header=False, compression='gzip')
    WHERE NOT starts_with(column0, '#')
  ", file_name)
  
  vcf_df <- DBI::dbGetQuery(con, query)
  
  # Clean up
  DBI::dbDisconnect(con)
  
  # Add column names to dataframe 
  # Read the file lines 
  lines <- readLines(gzfile(file_name), n = max_headlines)
  
  # Extract only lines starting with '#'
  header_lines <- grep("^#", lines, value = TRUE)
  
  # Get the last header line which contains column names (remove initial '#')
  header_line <- tail(header_lines, 1)
  
  # Remove the leading '#' and split by tab to get column names
  col_names <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  
  colnames(vcf_df) <- col_names
  
  # Return dataframe with vcf file info
  return(vcf_df)
}


#' @title compares the variants in a data.frame from a VCF file with variants in
#' another data.frame
#' @description Takes a data.frame with VCF variants and extracts those that are present
#' in another data.frame using DuckDB
#'
#' @param vcf_df a data.frame with variants from a VCF file
#' @param rsid_df a data.frame with genetic variants IDs 
#' @param memory_limit the number of GB set for DuckDB to operate (default: 2)
#' @param verbose if true, will print progress information (default: FALSE)
#' @return the "rsid_df" with merged columns containing information from the 
#' "gwas_df" 
#'
#' @examples
#' vargen_install("./vargen_data/")
#' DM1_simple <- vargen_pipeline(vargen_dir = "./vargen_data/", omim_morbid_ids = "222100",
#'                               fantom_corr = 0.25, outdir = "./", verbose = TRUE)
#' 
#' vcf_df <- load_vcf_file("./path_to_vcf_file,vcf.gz")
#' 
#' found_variants <- compare_vcf(vcf_df, DM1_simple)
#' @export
compare_vcf <- function(vcf_df, rsid_df, memory_limit = 2, verbose = FALSE){
  # Path to DuckDB file
  db_path <- "./duckdb_database.duckdb"
  
  # Open duckDB connection
  con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = db_path)
  
  # Set RAM usage limit to specified number of GB
  memory_query <- sprintf("SET memory_limit='%sGB'", memory_limit)
  DBI::dbExecute(con, memory_query)
  
  # Get dataframes
  duckdb::duckdb_register(con, "vcf_df", vcf_df)
  duckdb::duckdb_register(con, "rsid_df", rsid_df)
  
  # SQL query: join on df1.chrom = df2.chr and pos = pos
  query <- "
    SELECT df2.*, df1.REF, df1.ALT, df1.QUAL
    FROM vcf_df AS df1
    INNER JOIN rsid_df AS df2
    ON df2.chr = df1.CHROM AND df2.pos = df1.POS
  "
  
  # Run query and get result as dataframe
  intersect_df <- DBI::dbGetQuery(con, query)
  
  # Clean up
  DBI::dbDisconnect(con)
  
  # Return dataframe with shared variants in both dataframes
  return(intersect_df)
}