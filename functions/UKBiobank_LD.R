

LD.UKBiobank <- function(sumstats_path, 
                         results_path="./LD",
                         force_new_LD=F,
                         locus=NULL,
                         chimera=F, 
                         download_full_ld=F,
                         full_ld_path="./UKB_LD"){
   
  LD.UKB_find_ld_prefix <- function(chrom, min_pos){
    bp_starts <- seq(1,252000001, by = 1000000)
    bp_ends <- bp_starts+3000000
    i <- max(which(bp_starts<=min_pos)) 
    file.name <- paste0("chr",chrom,"_", bp_starts[i],"_", bp_ends[i])
    return(file.name)
  }
  
  
  LD.download_UKB_LD <- function(LD.file.list, 
                                 out.path="/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD/",
                                 alkes_url="https://data.broadinstitute.org/alkesgroup/UKBB_LD",
                                 background=T){
    flags <- ifelse(background,"-bqc","qc")
    dir.create(out.path, showWarnings = F, recursive = T)
    for(f in LD.file.list){ 
      gz.path <- file.path(alkes_url,paste0(f,".gz"))
      npz.path <- file.path(alkes_url,paste0(f,".npz"))
      # https://stackoverflow.com/questions/21365251/how-to-run-wget-in-background-for-an-unattended-download-of-files
      ## -bqc makes wget run in the background quietly
      system(paste("wget",gz.path,"-np",flags,"-P",out.path))
      cmd <- paste("wget",npz.path,"-np",flags,"-P",out.path)
      print(cmd)
      system(cmd)
    }
  }
  
  #### Support functions
  
  tryFunc <- function(input, func) {
    out <- tryCatch(
      {
        # Just to highlight: if you want to use more than one 
        # R expression in the "try" part then you'll have to 
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression 
        # in case the "try" part was completed successfully
        
        func(input) 
        # The return value of `readLines()` is the actual value 
        # that will be returned in case there is no condition 
        # (e.g. warning or error). 
        # You don't need to state the return value via `return()` as code 
        # in the "try" part is not wrapped insided a function (unlike that
        # for the condition handlers for warnings and error below)
      },
      error=function(cond) {
        message(paste("URL does not seem to exist:", input))
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      },
      warning=function(cond) {
        message(paste("URL caused a warning:", input))
        message("Here's the original warning message:")
        message(cond)
        # Choose a return value in case of warning
        return(NULL)
      },
      finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you 
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>' 
        message(paste("Processed URL:", input))
        message("Some other message at the end")
      }
    )    
    return(out)
  }
  
  
  printer <- function(..., v=T){if(v){print(paste(...))}}
  
  
  
  
  # Begin download
  finemap_DT <- data.table::fread(sumstats_path)
  alkes_url <- "https://data.broadinstitute.org/alkesgroup/UKBB_LD"
  URL <- alkes_url 
  chrom <- unique(finemap_DT$CHR)
  min_pos <- min(finemap_DT$POS) 
  file.name <- LD.UKB_find_ld_prefix(chrom=chrom, min_pos=min_pos)
  printer("+ UKB LD file name:",file.name)
  
  gz.path <- file.path(alkes_url, paste0(file.name,".gz")) 
  npz.path <- file.path(alkes_url, paste0(file.name,".npz"))
  
  chimera.path <- "/sc/orga/projects/pd-omics/tools/polyfun/UKBB_LD"
  UKBB.LD.file <- file.path(results_path, paste(locus,"UKB-LD.RDS",sep="_" ))
  
  if(download_full_ld){
    printer("+ LD:: Downloading full .gz/.npz UKB files and saving to disk.") 
    URL <- LD.download_UKB_LD(LD.file.list = file.name,
                               out.path = full_ld_path, 
                               background = F)
  } else {
    if(chimera){  
      if(file.exists(file.path(chimera.path, paste0(file.name,".gz")))  & 
         file.exists(file.path(chimera.path, paste0(file.name,".npz"))) ){
        printer("+ LD:: Pre-existing UKB LD gz/npz files detected. Importing...") 
        URL <- chimera.path
      }  
    } else {
      printer("+ LD:: Importing UKB LD file from alkesgroup database directly to R.") 
    }  
  }
  
  
  
  
  
  # RSIDs file
  # rsids <- data.table::fread(gz.path, nThread = 4) 
  if(file.exists(UKBB.LD.file) & force_new_LD!=T){
    printer("POLYFUN:: Pre-existing UKB_LD.RDS file detected. Importing...")
    LD_matrix <- readRDS(UKBB.LD.file)
  } else {
    # POLYFUN.load_conda(server = server)
    reticulate::source_python(file.path(polyfun,"load_ld.py"))
    printer("+ LD:: ...this could take some time...")
    
    ld.out <- tryFunc(input = file.path(URL, file.name), load_ld) 
    # LD matrix
    ld_R <- ld.out[[1]] 
    # head(ld_R)[1:10]
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
    printer("LD matrix dimensions", paste(dim(LD_matrix),collapse=" x "))
    printer("+ POLYFUN:: Saving LD =>",UKBB.LD.file)
    dir.create(dirname(UKBB.LD.file), showWarnings = F, recursive = T)
    saveRDS(LD_matrix, UKBB.LD.file)
  }  
  if(return_matrix){
    return(LD_matrix) 
  } else {
    return(UKBB.LD.file)
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

