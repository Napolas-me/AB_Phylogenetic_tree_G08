# Install Packages
install.packages("BiocManager")
BiocManager::install("msa")
BiocManager::install("ggtree")
BiocManager::install("ggplot2")
BiocManager::install("gridExtra")
install.packages("phangorn")
install.packages("stringdist")

# Load Packages
library("BiocManager")
library("Biostrings")
require("msa")
library("msa")

library("phangorn")
library("ape")
library("knitr")

library(ggtree)
library(ggplot2)
library(gridExtra)


# Tree Building Function
plot_tree = function(tree_plot, title_plot, max_x) {
  g = ggtree(tree_plot, color = "#00A499", size = 1)
  
  g = g + geom_tiplab(size = 4, color = "black", align = TRUE) +
    geom_nodepoint(size = 3, color = "#c7254e") +
    labs(title = title_plot, size = 6) +
    xlim(0, max_x) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = -5, b = 0, l = 0, unit = "pt"),
      legend.position = 'top'
    )
  
  g = g + geom_text(aes(x = branch, label = round(branch, 2)), size = 3,
                    color = "black", vjust = -0.5, hjust = 0)
  
  return(g)
}


### PREPARE AND PROCESS REQUIRED FILES ### 

# Path to the Folder with the FASTA Files **CHANGE TO DESIRED PATH**
folder_path <- "C:\\Mestrado\\1ano2semestre\\AlgoritmosBioinformatica\\Trabalho2\\Species\\Data"

# List All FASTA Files in the Folder (.fasta and .fa)
fasta_files <- list.files(path = folder_path, pattern = "\\.(fasta|fa)$", full.names = TRUE)

# Read FASTA Files into a DNAStringSet Object
sequences <- DNAStringSet()
for (file in fasta_files) {
  sequences <- c(sequences, readDNAStringSet(filepath = file))
}

# Convert Sequences to DNAStringSet Object
biostrings_sequences = DNAStringSet(sequences)



### MULTIPLE SEQUENCE ALLIGNMENT ###

# Start Multiple Sequence Alignment (Heads Up!!! Aprox 40 Minute Runtime)
msa_result_sample = msa(biostrings_sequences,method = "ClustalW")

# Obtain Sequence Names
sequence_names = rownames(msa_result_sample)

# Print Sequence Names, Sequence Content and MSA Results
print(sequence_names)
print(biostrings_sequences)
print(msa_result_sample)



### SAVE MSA RESULTS IN A FILE ###

# Convert the MSA Result to a String
msa_result_string <- as.character(msa_result_sample)

# Define the Output File Path **CHANGE TO DESIRED PATH**
output_file <- "C:\\Mestrado\\1ano2semestre\\AlgoritmosBioinformatica\\Trabalho2\\MSA.txt"

# Write the MSA Result to a Text File
write(msa_result_string, file = output_file)



### DISTANCE ALGORITHM ###

# Convert MSA Result to phyloDat Object
phyDat_msa_sample = as.phyDat(msa_result_sample)

# Add Custom Names (Using the Sequence Names)
names(phyDat_msa_sample) = c(sequence_names)

#Hamming Distance
D_hamming = dist.hamming(phyDat_msa_sample, ratio = FALSE)
D_hamming = round(as.matrix(D_hamming), 2)
print(D_hamming)



### TREE BUILDING ALGORITHMs ###

#Hamming Distance for Tree Algorithms
D_hamming = dist.hamming(phyDat_msa_sample)

# Nearest Neighbor Clustering (NNC)
tree_NNC = upgma(D_hamming, "single")
tree_NNC = midpoint(tree_NNC)

# Furthest Neighbor (FN)
tree_FN = upgma(D_hamming, "complete")
tree_FN = midpoint(tree_FN)

# Weighted Pair Group Method with Arithmetic Mean (WPGMA)
tree_WPGMA = wpgma(D_hamming)
tree_WPGMA = midpoint(tree_WPGMA)

# Unweighted Pair-Group Centroid Method (UPGMC)
tree_UPGMC = upgma(D_hamming, "centroid")
tree_UPGMC = midpoint(tree_UPGMC)

# Weighted Pair-Group Centroid Method (WPGGMC)
tree_WPGGMC = wpgma(D_hamming, "centroid")
tree_WPGGMC = midpoint(tree_WPGGMC)

# Unweighted Pair Group Method with Arithmetic Mean (UPGMA)
tree_UPGMA = upgma(D_hamming)
tree_UPGMA = midpoint(tree_UPGMA)

# Ward's Minimum Variance Method
tree_WARD = upgma(D_hamming, "ward")
tree_WARD = midpoint(tree_WARD)



### BUILDING THE PHYLOGENETIC TREES ###

# Nearest Neighbor Clustering (NNC) Tree
g_NNC = plot_tree(tree_NNC, "Nearest Neighbor Clustering (NNC)", 0.15)
grid.arrange(g_NNC)

# Furthest Neighbor (FN) Tree
g_FN = plot_tree(tree_FN, "Furthest Neighbor (FN)", 0.20)
grid.arrange(g_FN)

# Weighted Pair Group Method with Arithmetic Mean (WPGMA) Tree
g_WPGMA = plot_tree(tree_WPGMA, "Weighted Pair Group Method with Arithmetic Mean (WPGMA)", 0.17)
grid.arrange(g_WPGMA)

# Unweighted Pair-Group Centroid Method (UPGMC) Tree
g_UPGMC = plot_tree(tree_UPGMC, "Unweighted Pair-Group Centroid Method (UPGMC)", 0.10)
grid.arrange(g_UPGMC)

# Weighted Pair-Group Centroid Method (WPGGMC) Tree
g_WPGGMC = plot_tree(tree_WPGGMC, "Weighted Pair-Group Centroid Method (WPGGMC)", 0.10)
grid.arrange(g_WPGGMC)

# Unweighted Pair Group Method with Arithmetic Mean (UPGMA) Tree
g_UPGMA = plot_tree(tree_UPGMA, "Unweighted Pair Group Method with Arithmetic Mean (UPGMA)", 0.17)
grid.arrange(g_UPGMA)

# Ward's Minimum Variance Method Tree
g_ward = plot_tree(tree_WARD, "Ward's Minimum Variance Method", 0.6)
grid.arrange(g_ward)
