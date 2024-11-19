library(Biostrings)
library(dplyr)

# Step 1: Load the reference FASTA file
fasta_file <- "./input_data/reference_database.fasta"
fasta_sequences <- readAAStringSet(fasta_file)

fasta_df <- data.frame(
  ProteinID = names(fasta_sequences), 
  ProteinSequence = as.character(fasta_sequences),
  stringsAsFactors = FALSE
)
fasta_df <- fasta_df %>%
  mutate(ProteinID = sub(" .*", "", ProteinID)) 

# Step 2: Load the peptide dataset
peptide_df <- read.csv("./input_data/PTM_result.csv", stringsAsFactors = FALSE)

# Step 3: Extract naked peptide sequence and modification position
peptide_df <- peptide_df %>%
  mutate(
    NakedSequence = gsub("^[A-Z]\\.|\\.[A-Z]$", "", Peptides),    # Remove flanking residues
    NakedSequence = gsub(".-", "", NakedSequence),
    NakedSequence = gsub("-.", "", NakedSequence),
    NakedSequence = gsub("@", "", NakedSequence),
    ModificationIndices = gregexpr("&", NakedSequence), # Find position of '&'
    NakedSequence = gsub("&", "", NakedSequence)                 # Remove '&' for clean mapping
  )

# Step 4: Merge peptide data with FASTA protein sequences directly using ProteinID
merged_df <- peptide_df %>%
  left_join(fasta_df, by = c("Protein.accession" = "ProteinID"))

# Step 5: Map peptides to protein sequences and calculate modification sites
merged_df <- merged_df %>%
  rowwise() %>%
  mutate(
    PeptideStart = regexpr(NakedSequence, ProteinSequence)-1,  # Find peptide start position in protein
    ModificationSites = ifelse(
      PeptideStart > 0, 
      paste0(PeptideStart + unlist(ModificationIndices) - seq_along(unlist(ModificationIndices)), collapse = ";"),  # Adjust for '&' positions
      NA
    )
  ) %>%
  ungroup()

merged_df$ModificationSites <- as.character(merged_df$ModificationSites)
merged_df <- as.data.frame(merged_df)
merged_df <- merged_df[,c(1:3,7,8)]
writexl::write_xlsx(merged_df, "./output_data/modified_sites_in_proteins.xlsx")
