
source('./Rscripts/UKBiobank_LD.R')

ggLocusZoom <- function(sumstats_path, 
                        LD_path,
                        LD_units="r",
                        show_plot=T,
                        save_plot="./ggLocusZoom.png",
                        height=5, 
                        width=10,
                        title="",
                        leadSNP.line=T,
                        categorical_r2=T){
  library(ggbio) # This package ALIGNS lots of different data types in tracks (like UCSC Genome Browser).
  # See ggbio tutorial for more options: http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
  library(dplyr)
  library(EnsDb.Hsapiens.v75) # BiocManager::install("EnsDb.Hsapiens.v75")
  library(biovizBase)
  # library(emojifont)
  
  # ------------- [1.] IMPORT DATA ------------- 
  # Start with a dataframe containing your SNP-wise summary stats
  dat <- data.table::fread(sumstats_path)
  
  
  # ------------- [2. PREPARE DATA] ------------- 
  # Identify the lead SNP from the GWAS
  dat <- dat %>% mutate(leadSNP = ifelse(P==min(dat$P),T,F))
  
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
  } 
  # track.manhattan
  
  # [3.2] Manhattan track
  # You can pull gene models from EnsDb.Hsapiens.v75. 
  # But there's often way more than you can plot, so you have to filter them.
  # Get the row of the leadSNP so you can focus on gene annotations around it.
  lead.index <- which(gr.snp$leadSNP==T) 
  lw <- 10 
  track.genes <- autoplot(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, 
                          # Have to limit (can only handle depth < 1000)
                          which = gr.snp[(lead.index-lw):(lead.index+lw),], 
                          names.expr = "gene_name",  
                          aes(fill=gene_name, color=gene_name),
                          show.legend=T)  +  
    theme_bw()
   
  
  
  # [4.] ------------- MERGE TRACKS ------------- 
  # Add some features to tracks overall
  trks <- tracks(list("GWAS"=track.manhattan, "Genes"=track.genes), 
                 title = title, 
                 track.bg.color = "transparent",
                 track.plot.color = "transparent",
                 label.text.cex = .7, 
                 label.bg.fill = "grey12",
                 label.text.color = "white",
                 label.text.angle = 90) 
  #
  # BP_MB <- scales::trans_new("BP_MB", transform = function(x)x/1000000, inverse = function(x)x*1000000)
  trks <- trks +   
    scale_x_continuous(name = paste0("BP (",gr.snp$SEQnames[1],")"),
                       labels=function(x)x/1000000)
  if(leadSNP.line){
    trks <- trks +geom_vline(xintercept = lead.snp$POS, 
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
