#install.packages(c("data.table",
#                  "openxlsx",
#                   "ggplot2",
#                   "vegan",
#                   "zCompositions",
#                   "devtools",
#                   "propr",
#                   "partitions",
#                   "data.tree",
#                   "e1071"
#))

library(data.table) #work with dataframes
library(openxlsx) #open excel files


#################
#upload data
#################
path_to_counts_folder = "data"
counts_g = read.xlsx(file.path(path_to_counts_folder, "raw_counts_CRC.xlsx"),
                     sheet = "genus", rowNames = T) #first column from xlsx ro row names
metadata = fread(file.path(path_to_counts_folder, "metadata_CRC.txt"),
                 colClasses = "character")
alpha_diversity = read.xlsx(file.path(path_to_counts_folder, "alpha_diversity_CRC.xlsx"),
                            sheet = "shannon", rowNames = T)
counts_g = counts_g[metadata$sample,] #sorting counts_g based on sample column from metadata

# rename taxa to make names shorter
library(NearestBalance)
colnames(counts_g) <- make_nice_names(colnames(counts_g))

# filtering rare taxons:
# 3 or more counts in 40% of samples
include_taxon <- colSums(counts_g>2) > 0.4 * ncol(counts_g)
counts_filtered <- counts_g[, include_taxon]

# any sample with coverage < 3000? (Nope)
range(rowSums(counts_filtered))
coverage <- rowSums(counts_filtered)

# counting zeros
sum(counts_filtered == 0)/(ncol(counts_filtered)*nrow(counts_filtered))
sparsity(counts_filtered)


#################
# Create abundance table considering coverage and getting rid of zeros
#################
library(zCompositions)
abundance = cmultRepl(counts_filtered)
rowSums(abundance)




#configuring metadata
colnames(metadata)
names(metadata)[names(metadata) == "Diagnosis"] <- "case_control"
metadata_binary <- subset(metadata, case_control!="Large adenoma"&case_control!="Small adenoma")
abundance = abundance[metadata_binary$sample,]

# comparison of CRC-group with tumor-free control
library(NearestBalance)
nb_res <- nb_lm(abundance = abundance,
                metadata = metadata_binary,
                pred = "case_control",  # название столбца в метаданных
                type = "two_balances")

nb_res$nb 
nb_res$nb$impacts
nb_res$nb$b1
nb_res$lm_res
nb_res$sblm_summary
nb_res$noise #ideally <15% but we have 51%, more samples should be included

#manova for significance test
summary(manova(nb_res$lm_res)) #Pr(>F) <0.05 (yes, p = 0.000674)

#plot
require(balance)
library(ggplot2)
ilr <- balance.fromSBP(abundance, nb_res$nb$sbp)
ggplot(data.frame(ilr), aes(b1,b2, col = metadata_binary$case_control)) +
  geom_point() +
  coord_fixed() +
  theme_minimal()


#alpha diversity
alpha_diversity$sample <- rownames(alpha_diversity)
alpha <- merge(metadata, alpha_diversity, by.x = "sample", by.y = "sample")
alpha_crc <- subset(alpha, case_control=="Cancer")
alpha_normal <- subset(alpha, case_control=="Normal")

boxplot(alpha_normal$Shannon.index, alpha_crc$Shannon.index,
        main = "Alpha Diversity",
        names = c("Tumor-free", "Cancer")
)

#beta diversity
library(balance)
sbp <- sbp.fromRandom(abundance)
ilr <- balance.fromSBP(abundance, sbp)
d <- dist(ilr)

library(vegan)
perm <- how(nperm = 1000)
adonis2(d ~ metadata_binary$case_control, permutations = perm)



