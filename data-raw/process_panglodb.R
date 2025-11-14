library(orthogene)
library(dplyr)

pdb <- read.csv(
  file = "./data_raw/PanglaoDB_markers_27_Mar_2020.tsv",
  sep = "\t",
  header = TRUE
)

colnames(pdb) <- gsub("\\.", "_", colnames(pdb))

panglaodb_hs <- pdb[grepl("Hs", pdb$species), ] %>%
  rename(gene = "official_gene_symbol") %>%
  rename(cell_type = "cell_type")

panglaodb_mm <- pdb[grepl("Mm", pdb$species), ] %>%
  rename(gene = "official_gene_symbol") %>%
  rename(cell_type = "cell_type")

gene_df <- orthogene::convert_orthologs(
  gene_df = panglaodb_mm,
  gene_input = "gene",
  gene_output = "columns",
  input_species = "human",
  output_species = "mouse",
  non121_strategy = "drop_both_species",
  method = "gprofiler"
)

panglaodb_mm <- panglaodb_mm[which(panglaodb_mm$gene %in% gene_df$input_gene), ]
if (identical(panglaodb_mm$gene, gene_df$input_gene)) {
  gene_df <- gene_df[match(panglaodb_mm$gene, gene_df$input_gene), ]
}
panglaodb_mm$gene <- gene_df$ortholog_gene

usethis::use_data(panglaodb_hs, overwrite = TRUE)
usethis::use_data(panglaodb_mm, overwrite = TRUE)
