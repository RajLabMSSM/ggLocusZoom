

LD.UKBiobank <- function(sumstats_path=NULL,
                         chrom=NULL,
                         min_pos=NULL,
                         output.path="./LD",
                         force_new_LD=F,
                         locus=NULL,
                         chimera=F, 
                         download_full_ld=F, 
                         method="axel", 
                         nThreads=4,
                         return_matrix=F,
                         remove_tmps=F){
  # Find LD prefix
  LD.UKB_find_ld_prefix <- function(chrom, min_pos){
    bp_starts <- seq(1,252000001, by = 1000000)
    bp_ends <- bp_starts+3000000
    i <- max(which(bp_starts<=min_pos)) 
    file.name <- paste0("chr",chrom,"_", bp_starts[i],"_", bp_ends[i])
    return(file.name)
  }
   
  
  LD.download_UKB_LD <- function(LD.prefixes, 
                                 output.path="/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD/",
                                 alkes_url="https://data.broadinstitute.org/alkesgroup/UKBB_LD",
                                 background=T, 
                                 force_overwrite=F,
                                 method="axel"){ 
    for(f in LD.prefixes){ 
      gz.url <- file.path(alkes_url,paste0(f,".gz"))
      npz.url <- file.path(alkes_url,paste0(f,".npz")) 
       
      for(furl in c(gz.url, npz.url)){
        if(tolower(method)=="axel"){
          out.file <- axel(input.url = furl,
                           output.path = output.path,
                           background = background,
                           nThreads = nThreads, 
                           force_overwrite = force_overwrite)
        }
        if(tolower(method)=="wget"){
          out.file <- wget(input.url = furl, 
                           output.path = output.path, 
                           background = background,
                           force_overwrite = T) 
        } 
      }
    }
    return(gsub("*.npz$","",out.file))
  }
  
  #### Support functions
  
  tryFunc <- function(input, func, ...){
    out <- tryCatch(
      { 
        func(input, ...)  
      },
      error=function(cond) {
        message(paste("URL does not seem to exist:", input))
        message("Here's the original error message:")
        message(cond) 
        return(NA)
      },
      warning=function(cond) {
        message(paste("URL caused a warning:", input))
        message("Here's the original warning message:")
        message(cond) 
        return(NULL)
      },
      finally={ 
        message(paste("Processed URL:", input))
        message("Some other message at the end")
      }
    )    
    return(out)
  }
  
  printer <- function(..., v=T){if(v){print(paste(...))}}
  
  ## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ## Quickstart
  # sumstats_path="./example_data/BST1_Nalls23andMe_2019_subset.txt"
  # chrom=NULL
  # min_pos=NULL
  # output.path="./LD"
  # force_new_LD=F
  # locus="BST1"
  # chimera=F
  # download_full_ld=T
  # method="axel"
  # nThreads=4
  
  # Begin download
  if(!is.null(sumstats_path)){
    printer("+ Assigning chrom and min_pos based on summary stats file")
    # sumstats_path="./example_data/BST1_Nalls23andMe_2019_subset.txt"
    finemap_DT <- data.table::fread(sumstats_path, nThread = 4)
    chrom <- unique(finemap_DT$CHR)
    min_pos <- min(finemap_DT$POS) 
    LD.prefixes <- LD.UKB_find_ld_prefix(chrom=chrom, min_pos=min_pos)
  }
  
  alkes_url <- "https://data.broadinstitute.org/alkesgroup/UKBB_LD"
  URL <- alkes_url 
  
  printer("+ UKB LD file name:",LD.prefixes)  
  chimera.path <- "/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD"
  UKBB.LD.RDS <- file.path(output.path, paste(locus,"UKB-LD.RDS",sep="_" ))
  
  if(download_full_ld){
    printer("+ LD:: Downloading full .gz/.npz UKB files and saving to disk.") 
    URL <- LD.download_UKB_LD(LD.prefixes = LD.prefixes,
                              output.path = output.path, 
                              background = F,
                              force_overwrite = force_new_LD, 
                              method = method)
    server <- F
  } else {
    if(chimera){  
      if(file.exists(file.path(chimera.path, paste0(LD.prefixes,".gz")))  & 
         file.exists(file.path(chimera.path, paste0(LD.prefixes,".npz"))) ){
        printer("+ LD:: Pre-existing UKB LD gz/npz files detected. Importing...") 
        URL <- chimera.path
      }  
    } else {
      printer("+ LD:: Importing UKB LD file from alkesgroup database directly to R.") 
    }  
  }
  
  
  # RSIDs file
  # rsids <- data.table::fread(gz.path, nThread = 4) 
  if(file.exists(UKBB.LD.RDS) & force_new_LD!=T){
    printer("+ LD:: Pre-existing UKB_LD.RDS file detected. Importing...")
    LD_matrix <- readRDS(UKBB.LD.RDS)
  } else { 
    printer("+ LD:: ...this could take some time...")
    reticulate::source_python(file.path("functions","load_ld.py"))
    # load_ld(ld_prefix=URL, server=F)
    ld.out <- tryFunc(input = URL, load_ld, server = server)
    # LD matrix
    ld_R <- ld.out[[1]] 
    # head(ld_R)[1:10,]
    # SNP info 
    # ld_snps <- data.table::data.table( reticulate::py_to_r(ld.out[[2]]) ) 
    ld_snps <- data.table::data.table(ld.out[[2]])
    row.names(ld_R) <- ld_snps$rsid
    colnames(ld_R) <- ld_snps$rsid 
    
    # remove(ld.out) 
    # ld_snps.sub <- subset(ld_snps, position %in% finemap_DT$POS)
    indices <- which(ld_snps$position %in% finemap_DT$POS)
    ld_snps.sub <- ld_snps[indices,]
    LD_matrix <- ld_R[indices, indices]
    row.names(LD_matrix) <- ld_snps.sub$rsid
    colnames(LD_matrix) <- ld_snps.sub$rsid 
    LD_matrix[is.na(LD_matrix)] <- 0
    # Save LD matrix as RDS
    printer("LD matrix dimensions", paste(dim(LD_matrix),collapse=" x "))
    printer("+ POLYFUN:: Saving LD =>",UKBB.LD.RDS)
    dir.create(dirname(UKBB.LD.RDS), showWarnings = F, recursive = T)
    saveRDS(LD_matrix, UKBB.LD.RDS)
    
    if(remove_tmps){
      printer("+ Removing .gz/.npz files.")
      file.remove(paste0(URL,".gz"))
      file.remove(paste0(URL,".npz"))
    }
  }  
  if(return_matrix){
    return(LD_matrix) 
  } else {
    return(UKBB.LD.RDS)
  } 
}





# LD.UKBiobank_multi_download <- function(out.path = "/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD/"){
#   # Download all UKBB LD files
#   # wget -r -np -A '*.gz' https://data.broadinstitute.org/alkesgroup/UKBB_LD/
#   # wget -r -np -A '*.npz' https://data.broadinstitute.org/alkesgroup/UKBB_LD/
#   
#   # Figure out the names of LD files you'll need
#   merged_DT <- merge_finemapping_results(minimum_support = 0)
#   locus_coords <- merged_DT %>% dplyr::group_by(Gene) %>% 
#     dplyr::summarise(chrom=unique(CHR), min_pos=min(POS), max_pos=max(POS))
#   LD.file.list <- lapply(1:nrow(locus_coords), function(i){
#     return(LD.find_ld_prefix(chrom = locus_coords$chrom[i], min_pos=locus_coords$min_pos[i]))
#   }) %>% unlist()
#   
#   LD.download_UKB_LD(LD.file.list = LD.file.list,
#                      out.path = out.path)
#   return(LD.file.list)
# }


# POLYFUN.load_conda <- function(server=F){
#   printer("POLYFUN:: Activating polyfun_venv...")
#   if(server){ 
#     reticulate::use_condaenv("polyfun_venv")
#   } else {
#     reticulate::use_condaenv("polyfun_venv", conda = "/usr/local/anaconda3/condabin/conda")
#     # conda_list("/usr/local/anaconda3/condabin/conda")
#   } 
# }
# 


