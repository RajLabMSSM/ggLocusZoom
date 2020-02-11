
source('./functions/UKBiobank_LD.R')
source('./functions/downloaders.R')

# BiocManager::install("interactiveDisplay")
# https://www.bioconductor.org/packages/release/bioc/vignettes/interactiveDisplay/inst/doc/interactiveDisplay.pdf

# Support functions for converting gene names
hgnc_to_ensembl <- function(gene_symbols){  
  # columns(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                   keys = gene_symbols,
                                   keytype = "SYMBOL",
                                   column = "GENEID")
  return(conversion)
}
 
ensembl_to_hgnc <- function(ensembl_ids){
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                   keys = ensembl_ids,
                                   keytype = "GENEID",
                                   column = "SYMBOL")
  return(conversion)
}



ggLocusZoom <- function(sumstats_path, 
                        LD_path,
                        LD_units="r",
                        show_plot=T,
                        save_plot="./ggLocusZoom.png",
                        height=5, 
                        width=10,
                        title="",
                        leadSNP.line=T,
                        categorical_r2=T,
                        index_snp="leadGWAS",
                        xlims=NULL, 
                        max_transcripts=3){
  library(ggbio) # This package ALIGNS lots of different data types in tracks (like UCSC Genome Browser).
  # See ggbio tutorial for more options: http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
  library(dplyr)
  library(EnsDb.Hsapiens.v75) # BiocManager::install("EnsDb.Hsapiens.v75")
  library(biovizBase)
  # library(emojifont)
  
  # ------------- [1.] IMPORT DATA ------------- 
  # Start with a dataframe containing your SNP-wise summary stats
  dat <- data.table::fread(sumstats_path, nThread = 4)
  
  
  # ------------- [2. PREPARE DATA] ------------- 
  # Identify the lead SNP from the GWAS
  if(index_snp=="leadGWAS"){
    dat <- dat %>% mutate(leadSNP = ifelse(P==min(dat$P),T,F))
  } else {
    dat <- dat %>% mutate(leadSNP = ifelse(SNP==index_snp,T,F))
  }
 
  
  # Get your LD 
  if(endsWith(LD_path,".RDS")){
    LD_matrix <- readRDS(LD_path)
  } else {
    LD_matrix <- data.table::fread(LD_path)
  }
 
  ## Get the r2 with your lead SNP
  r2.dat <- LD_matrix[,subset(dat,leadSNP==T)$SNP] 
  if(LD_units=="r") r2.dat <- r2.dat^2
  r2.dat <- data.frame(SNP=names(r2.dat), r2=r2.dat)
  ### Turn r2 into categories
  r2.dat <- r2.dat %>% mutate(r2.group=cut(r2, breaks =  c(0,.2,.4,.6,.8,1), right=F, 
                                 labels=c("very low","low","medium","high","very high")))
  ## Merge the r2 with your original
  dat <- merge(dat, r2.dat, by="SNP", all.x=T)  
  
  # Now, convert it to a GRanges object
  dat <- dat %>% dplyr::mutate(SEQnames = paste0("chr",CHR))
  gr.snp <- biovizBase::transformDfToGr(dat, seqnames = "CHR", start = "POS", end = "POS")
  
  lead.snp <- data.frame(subset(gr.snp, leadSNP))
  
  # ------------- [3.] CREATE TRACKS ------------- 
  
  # [3.1] Manhattan track 
  ## Replace r2.group with r2 for continous LD
  if(categorical_r2){
    track.manhattan <- plotGrandLinear(obj=gr.snp, aes(y=-log10(P), color=r2.group, shape=leadSNP),
                                       show.legend=T, alpha=.75) + 
      scale_color_brewer(palette = "Spectral",direction = -1) + 
      scale_shape_manual(values=c(16,18)) +
      annotate(geom="text", x = lead.snp$POS, y = -log10(lead.snp$P)*1.05,label=lead.snp$SNP, size=3) +
      theme_bw()
  } else {
    track.manhattan <- plotGrandLinear(obj=gr.snp, aes(y=-log10(P), color=r2, shape=leadSNP),
                                       show.legend=T, alpha=.75) + 
      scale_color_gradient(low="blue", high="red", limits = c(0,1)) + 
      scale_shape_manual(values=c(16,18)) +
      annotate(geom="text", x = lead.snp$POS, y = -log10(lead.snp$P)*1.05,label=lead.snp$SNP, size=3) +
      theme_bw()
    # ggbio::zoom_in(fac = 10) 
  } 
  # track.manhattan
  
  # [3.2] Manhattan track
  # You can pull gene models from EnsDb.Hsapiens.v75. 
  # But there's often way more than you can plot, so you have to filter them.
  # Get the row of the leadSNP so you can focus on gene annotations around it.
  lead.index <- which(gr.snp$leadSNP==T)  
  # if(gene_level==T){
  #   print("+ Annotating at gene-level.")
  #   db <- ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  #   db.gr <- data.table::as.data.table(db) %>%
  #     dplyr::mutate(index=row.names(.))  
  #   db.gr <- db.gr %>% dplyr::group_by(gene_name) %>% dplyr::slice(1)
  # }  
    print("+ Annotating at transcript-level.") 
    db.gr <- ensembldb::transcripts(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) %>%
      data.table::as.data.table() %>%
      dplyr::mutate(index=row.names(.)) %>% 
      dplyr::group_by(gene_id) %>% 
      dplyr::slice(1:max_transcripts)  
  if(!"symbol" %in% colnames(db.gr)){
    db.gr$symbol <- ensembl_to_hgnc(db.gr$gene_id)
  }
  # Subset to only the region encompassed by the sumstats
  db.gr <- subset(db.gr, seqnames == unique(dat$CHR) &
                    start >= min(dat$POS) &
                    end <= max(dat$POS)) %>%  
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) 
  db.gr$symbol <- factor(db.gr$symbol, levels = unique(db.gr$symbol), ordered = T)
    # dplyr::group_by(gene_name) %>%
    # dplyr::top_n(n = top_transcripts, wt=seqnames)   
  edb <- addFilter( EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,  AnnotationFilter::TxIdFilter(db.gr$tx_id))   
  track.genes <- autoplot(edb, 
                          # Have to limit (can only handle depth < 1000)
                          which = db.gr,
                          names.expr = "gene_name",  
                          aes(fill=gene_name, color=gene_name),
                          show.legend=T)  +  
    theme_bw()
   
  
  if(is.null(xlims)){
    xlims <- c(min(dat$POS), max(dat$POS))
  }
  # [4.] ------------- MERGE TRACKS ------------- 
  # Add some features to tracks overall
  trks <- tracks(list("GWAS"=track.manhattan, "Genes"=track.genes), 
                 title = title, 
                 track.bg.color = "transparent",
                 track.plot.color = "transparent",
                 label.text.cex = .7, 
                 label.bg.fill = "grey12",
                 label.text.color = "white",
                 label.text.angle = 90,
                 xlim = xlims) 
  #
  # BP_MB <- scales::trans_new("BP_MB", transform = function(x)x/1000000, inverse = function(x)x*1000000)
  trks <- trks +   
    scale_x_continuous(name = paste0("BP (",gr.snp$SEQnames[1],")"),
                       labels=function(x)x/1000000)
  if(leadSNP.line){
    trks <- trks + geom_vline(xintercept = lead.snp$POS, 
                              color="red", alpha=.6, size=.3, linetype='solid') 
  }
  # Print
  if(show_plot) print(trks)
  # save
  if(save_plot!=F){
    print(paste("Saving plot ===>",save_plot))
    ggsave(filename = save_plot, plot =  trks, dpi=300, height=height, width=width)
  }
  return(trks)
}
