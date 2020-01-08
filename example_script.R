
source("./functions/ggLocusZoom.R")


LD_path <- LD.UKBiobank(# Specify where the summary stats file is.
  sumstats_path="./example_data/BST1_Nalls23andMe_2019_subset.txt", 
  
  # The folder where you want to save the LD matrix 
  # (defaults to ./LD)
  results_path="./LD",
  
  # If =T , will overwrite a pre-existing LD file with the same name.
  force_new_LD=F,
  
  # The same of the locus 
  # (defaults to "_")
  locus="BST",
  
  # Use pre-downloaded LD files on Chimera computing cluster 
  # (Mount Sinai employees and affiliates only). 
  # Must be logged onto Chimera and have 
  # access to the 'pd-omics' project.
  chimera=F, 
  
  # [** WARNING **]: Only change these defaults if you have plenty of extra storage. Each of these files is ~1-3GB.
  ## Download and save the full .npz/.gz files.
  download_full_ld=F,
  ## Specify where to save these files.
  full_ld_path="./UKB_LD")
print(LD_path)






gglz <- ggLocusZoom(# Specify where the summary stats file is.
  sumstats_path="./example_data/BST1_Nalls23andMe_2019_subset.txt",
  
  # Specify where the pre-computed LD matrix is.
  LD_path="../example_data/BST1_UKB-LD.RDS",
  
  # Is the LD matrix in units of r? (as opposed to r2)
  LD_units="r",
  
  # Show the plot when it's ready?
  show_plot=T,
  
  # Save the plot? (=F if you don't want to)
  save_plot="./BST1_ggLocusZoom.png",
  
  # Saved plot height.
  height=5, 
  
  # Saved plot width.
  width=10,
  
  # Optional plot title.
  title="",
  
  # Draw a vertical line at the position of the lead SNP.
  leadSNP.line=T,
  
  # Plot LD with the leda SNP as a categorical variable 
  # (very low, low, medium, high, very high) or a continuous variable (0-1).
  categorical_r2=T)