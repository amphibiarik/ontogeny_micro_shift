## Code used to perform analyses and generate figures comparing bacterial diversity of the salamander microbiome through development from Hartmann et al. 2023, in *Microbiology*
### Publication is open access available at
https://doi.org/10.1099/mic.0.001399

# This repository contains:

| Script  | Description |
| ------------- | ------------- |
| Hartmann_et_al_Ontogeny.R  | Set up, load libraries, phyloseq creation, analyses, and figures  |

| File  | Description |
| ------------- | ------------- |	
| MappingFileSN.csv  | Sample data for newts used in microbiome study  |
| MappingFileSN.txt  | Sample data for newts used in microbiome study, used to map samples to phyloseq object  |
| Seq.tab.nochim.rds | Sequence table to with chimeric reads removed. Used to create phyloseq object  |
| Taxa_taxid.rds | Microbial taxa table generated using SILVA database. Used to assign bacterial taxa to phyloseq object  |
| treeNJ.rds | Phylogenetic tree of bacterial taxa used for phyloseq generation  |
| res_W_ANCOM.LS.csv | Results of ANCOMBC by newt life stage |
