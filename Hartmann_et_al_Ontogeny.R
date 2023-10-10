##### Code to reproduce analyses from Hartmann et al. Microbiology
#Title: Ontogeny drives shifts in skin bacterial communities in facultatively paedomorphic salamanders
#Please contact corresponding author Arik Hartmann for questions regarding the code and data (arikhartmann.ufl.edu)


library("devtools")
library("dada2")
library("DECIPHER")
library(phyloseq)
library(phangorn)
library(ggplot2)
library("BiocGenerics")
library(reshape2)
library(patchwork)
library(microbiome)
library(qiime2R)
library(vegan)
library(MetBrewer)
library(microbiomeutilities)
library(ggpubr)
library(gridExtra)
library(picante)
library(ANCOMBC)
library(MicEco)
library(ecole)

## Summarized sample information
sample_data <- read.csv("MappingFileSN.csv")


#### running with phylotree ####
treeNJ <- readRDS(file='treeNJ.rds')

seqtab.nochim <- readRDS("seqtab.nochim.rds")

# add taxa
taxa <- readRDS("taxa_taxid.rds")

#generate phyloseq object with phylotree
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa), phy_tree(treeNJ))

# add sample data
map <- import_qiime_sample_data("MappingFileSN.txt")
rownames(map) <- map$sample.id

#merge phyloseq with map
ps <- merge_phyloseq(ps,map)
ps
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# set seed and root tree
set.seed(711)

phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))

#Remove controls and prune low reads 
ps1 <- prune_samples(sample_names(ps) != "43-SNZymo-std1-515rcbc29", ps) # control
ps1 <- prune_samples(sample_names(ps1) != "46-SNNA2-515rcbc65", ps1) # low reads
ps1 <- prune_samples(sample_names(ps1) != "42-SN1066-515rcbc17", ps1) # low reads
# this sample was dropped in the dirt :(
ps1 <- prune_samples(sample_names(ps1) != "37-SN1051-515rcbc52", ps1) 


# Select only Bacteria that can be identified to the Phyla level
bacteria <- subset_taxa(ps1, domain=="Bacteria" & phylum!="Unassigned")
bacteria

# identify lowest reads of retained samples to set rarification depth
summary(sample_sums(bacteria))

bact.rar <- rarefy_even_depth(bacteria, sample.size = 37897)

# Get metadata from PS object
bac.meta <- meta(bact.rar)

# Add the rownames as a new colum for easy integration later.
bac.meta$sam_name <- rownames(bac.meta)

## Alpha diversity metrics ####

# calculate Alpha Diversity measures
bac.div <- alpha(bact.rar, index = "all")
# Add the rownames to diversity table
bac.div$sam_name <- rownames(bac.div)

# for Faith's PD
bact.otu <- as.data.frame(bact.rar@otu_table)
bact.tree <- bact.rar@phy_tree

bact.pd <- pd(bact.otu, bact.tree, include.root=T) 

bac.meta$Phylogenetic_Diversity <- bact.pd$PD

#set order for infetcion and life stages
bac.meta$Infection_type = factor(bac.meta$Infection_type, levels = c("Bd", "Coinfected", "Rv", "None"))
bac.meta$life_stage = factor(bac.meta$life_stage, levels = c("L", "PA", "MA", "A"))

# merge bacterial diversity with meta data
div.df <- merge(bac.div,bac.meta, by = "sam_name")

# check the tables
colnames(bac.meta)
colnames(div.df)


# make a pairwise list that we want to compare.
inf.pairs <- list(c("Bd", "Rv"), c("Bd", "Coinfected"), c("Bd", "None"), c("Rv", "Coinfected"), c("Rv", "None"),
                       c("Coinfected", "None"))

ls.pairs <- list(c("PA", "L"), c("PA", "MA"), c("PA", "A"), c("MA", "L"), c("MA", "A"),
                c("L", "A"))

# set significant cutoff
pval <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("****", "***", "**", "*", "n.s")
)

### GLMs for alpha diversity####

#Faith's PD
glmpd <-glm(Phylogenetic_Diversity~life_stage+Infection_type+Pond+SVL_cm+Mass+Date, data=div.df, family=gaussian)
summary(glmpd)
anova.glm(glmpd, test="F")

#Shannon
glmshan <-glm(diversity_shannon~Infection_type+life_stage+SVL_cm+Mass+Date+Pond, data=div.df)
summary(glmshan)
anova.glm(glmshan, test="F")

# richness
glmrich <- glm(log(observed)~life_stage+Infection_type+SVL_cm+Mass+Date, data=div.df, family=gaussian)
shapiro.test(rstandard(lm_s1))
summary(glmrich)
anova.glm(glmrich, test="F")

##Beta Diversity pairwise PERMANOVAs ####
# Calculate Bray-curtis, unweighted, and weighted UNIfrac distances
bray_ls <- phyloseq::distance(bact.rar, "bray")
uw_uf <- phyloseq::distance(bact.rar, "unifrac")
w_uf <-phyloseq::distance(bact.rar, "wunifrac")

#Life stage - Pairwise PERMANOVA followed by HSD 
#Bray-curtis distance
bray_pw <- permanova_pairwise(x=bray_ls, grp=div.df$life_stage, permutations=999)
bray_pw
bc.disper <- betadisper(bray_ls, div.df$life_stage)
anova(bc.disper)
TukeyHSD(bc.disper)

#Unweighted UniFrac
uwuf_pw <- permanova_pairwise(x=uw_uf, grp=div.df$life_stage, permutations = 999)
uwuf_pw
uwuf.disper <- betadisper(uw_uf, div.df$life_stage)
anova(uwuf.disper)
TukeyHSD(uwuf.disper)

#Weighted UniFrac
wuf_pw <- permanova_pairwise(x=w_uf, grp=div.df$life_stage, permutations = 999)
wuf.disper <- betadisper(w_uf, div.df$life_stage)
anova(wuf.disper)
TukeyHSD(wuf.disper)

#Infection - Pairwise PERMANOVA, PERMDISP + HSD

#Bray-curtis distance
bray_inf <- permanova_pairwise(x=bray_ls, grp=div.df$Infection_type, permutations=999)
bray_inf
bc.inf.disper <- betadisper(bray_ls, div.df$Infection_type)
anova(bc.inf.disper)
TukeyHSD(bc.inf.disper)

#Unweighted UniFrac
uwuf_inf <- permanova_pairwise(x=uw_uf, grp=div.df$Infection_type, permutations=999)
uwuf_inf
uwuf.inf.disper <- betadisper(uw_uf, div.df$Infection_type)
anova(uwuf.inf.disper)
TukeyHSD(uwuf.inf.disper)

#Weighted UniFrac
wuf_inf <- permanova_pairwise(x=w_uf, grp=div.df$Infection_type, permutations=999)
wuf_inf
wuf.inf.disper <- betadisper(w_uf, div.df$Infection_type)
anova(wuf.inf.disper)
TukeyHSD(wuf.inf.disper)

#### Figure 1 ####
#Set palette for life stages
lspal <- met.brewer("Greek", n=5)[c("L"=1,"PA"=2,"MA"=3,"A"=5)]

#Fig 1A - Faith's Phylogenetic Diversity
#
pd.plot <- ggboxplot(bac.meta,
                     x = "life_stage",
                     y = "Phylogenetic_Diversity",
                     fill = "life_stage",
                     palette = "jco",
                     ylab = "Phylogenetic Diversity",
                     xlab = "Life stage",
                     legend = "right", alpha=0.75
)
pd.plot <- pd.plot + rotate_x_text()

faith_plot <- pd.plot 

faith <- faith_plot + scale_fill_manual(name="Life stage", values=lspal,
                                        labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"))+
  labs(x=NULL, y="Faith's PD\n")+
  scale_x_discrete(labels = c('Larva','Paedomorph','Metamorphosing', "Adult"))+
  ylim(0, 70)+
  stat_compare_means(
    comparisons = inf.pairs, bracket.size = 0.9, size=8,
    v.just=0.25,
    label = "p.signif",
    label.y = c(NA, 48, 54, NA, NA, 60))+
  theme_classic() + theme(axis.line = element_line(size=0.75),
                          axis.ticks.length = unit(0.25, "cm"),
                          axis.ticks = element_line(linewidth = 0.75),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_text(size = 20, color="black"),
                          axis.text.y = element_text(size=20, color="black"),
                          axis.title = element_text(size=24),
                          legend.title = element_blank(),
                          legend.position = "none")
plot(faith)

#Fig 1B - Bray-Curtis PCo1

# transform Phyloseq object to gain relative abundance
bact.rar.rel <- microbiome::transform(bact.rar, "compositional")

# Ordinate relative abundances by Bray-curtis distance
pcoa_bc = ordinate(bact.rar.rel, "PCoA", "bray") 

bc <- plot_ordination(bact.rar.rel, pcoa_bc, color="life_stage") + 
  geom_point(aes(fill=life_stage, color=life_stage), size = 5, alpha =0.75)+
  scale_fill_manual(name="Life stage", values=lspal, labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"))+
  scale_color_manual(name="Life stage", values=lspal, labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"))
  
bc <- bc + stat_ellipse(aes(group=life_stage, fill=factor(life_stage)), 
                        geom="polygon", alpha=0.2, show.legend = F, linetype="blank")
# remove initial ordination layer to eliminate double dots when using a transparent ggplot
bc$layers <- bc$layers[-1]
bc

# add aesthetics
bc_final <- bc + theme_classic() + theme(axis.line = element_line(size=0.75),
                                         axis.ticks.length = unit(0.25, "cm"),
                                         axis.ticks = element_line(linewidth = 0.75),
                                         axis.text = element_text(size = 22, color="black"),
                                         axis.title = element_text(size=24),
                                         legend.position = "none")

bc_final

# repeat for Richness
rich.plot <- ggboxplot(div.df,
                       x = "life_stage",
                       y = "observed",
                       fill = "life_stage",
                       palette = "jco",
                       ylab = "ASV richness",
                       xlab = "Life stage",
                       legend = "right", alpha=0.75
)
rich.plot <- rich.plot + rotate_x_text()

rich <- rich.plot + scale_fill_manual(name="Life stage", values=lspal,
                                      labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"))+
  labs(x=NULL, y="ASV richness\n")+
  scale_x_discrete(labels = c('Larva','Paedomorph','Metamorphosing', "Adult"))+
  ylim(0, 900)+
  stat_compare_means(
    comparisons = ls.pairs,
    bracket.size = 0.9, size=8,
    v.just=0.25,
    label = "p.signif",
    label.y = c(NA, 740, 800, NA, NA, 860))+
  theme_classic() + theme(axis.line = element_line(size=0.75),
                          axis.ticks.length = unit(0.25, "cm"),
                          axis.ticks = element_line(linewidth = 0.75),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_text(size = 20, color="black"),
                          axis.text.y = element_text(size=20, color="black"),
                          axis.title = element_text(size=24),
                          legend.title = element_text(size=20),
                          legend.text = element_text(size=18))
plot(rich)

#### Figure 2 - Core Bacteria #### 

#Fig 2A - Core All
core_all <- ps_venn(bact.rar,
                    group="life_stage",
                    quantities=list(type=c("counts","percent"), font=3, round=2, cex=1),
                    fraction = 0.00,
                    weight =  F,
                    relative = T,
                    plot = T,
                    main=list(labels=c("Core-All"), font=1, cex=3),
                    labels=list(labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"), font=1, cex=1.5),
                    fill=lspal, alpha=0.40)

#Fig 2B - Core75
core75 <- ps_venn(bact.rar,
                  group="life_stage",
                  quantities=list(type=c("counts","percent"),font=3, round=2, cex=1),
                  fraction = 0.75,
                  weight =  F,
                  relative = T,
                  plot = T,
                  main=list(labels=c("Core75"), font=1, cex=3),
                  labels=list(labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"),font=1, cex=1.5),
                  fill=lspal, alpha=0.60)

#Fig 2C - Core90
core90 <- ps_venn(bact.rar,
                  group="life_stage",
                  quantities=list(type=c("counts","percent"),font=3, round=2, cex=1),
                  fraction = 0.90,
                  weight =  F,
                  relative = T,
                  plot = T,
                  main=list(labels=c("Core90"), font=1, cex=3),
                  labels=list(labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"),font=1, cex=1.5),
                  fill=lspal, alpha=0.75)

# Arrange plot layout and annotate each figure
Fig2 <- (wrap_elements(core_all) + wrap_elements(core75) + wrap_elements(core90)) +
  plot_layout(ncol=3)+
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size=22))

Fig2

ggsave("Fig2.tiff", plot=Fig2, units = "in", width = 21, height = 8, dp = 300,compression = 'lzw')

####Figure 3 - Relative abundance of bacterial families ####

# set palette with >20 colors
fpal20 <- met.brewer("Renoir", n=21)

Fam20 <- plot_composition(Fam20,
                           average_by = "life_stage", 
                           transform = "compositional",
                           fill="family")

Fam20 <- Fam20 + scale_fill_manual(name="Bacterial family", values=fpal20)


Fam20 <- Fam20 + labs(x=NULL, y="Relative abundance\n", title=NULL)+theme_classic()+
  scale_x_discrete(labels=c("Larva", "Paedomorph", "Metamorphosing", "Adult"))+
  theme(axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=22, color="black"),
        axis.title.y = element_text(size=24),
        axis.line = element_line(size=0.75),
        axis.ticks.y = element_line(size=0.75),
        axis.ticks.length.y = unit(0.25, "cm"),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=20, color="black"),
        legend.text = element_text(size=18))+guides(fill=guide_legend(ncol=1))

ggsave("Fig3.tiff", pcomp2, dpi=300, height = 6, width = 9, units='in')

#### ANCOMBC ####

# ANCOM at family level
pseq_fam = aggregate_taxa(bact.rar, "family")
tax_mat = as(tax_table(pseq_fam), "matrix")

out = ancombc(
  phyloseq = pseq_fam, 
  formula = "life_stage", 
  p_adj_method = "fdr", 
  zero_cut = 0.90, # by default prevalence filter of 10% is applied
  lib_cut = 0, 
  group = "life_stage", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
resfam <- out$res

res.global.fam <- out$res_global

#separate for each life stage, rename and combine
write.csv(resfam$se, "res_se_ANCOM.csv", row.names = T)
write.csv(resfam$diff_abn, "res_diff_abn_ANCOM.csv", row.names = T)
write.csv(resfam$p_val, "res_p_ANCOM.csv", row.names = T)
write.csv(resfam$q_val, "res_q_ANCOM.csv", row.names = T)
write.csv(resfam$W, "res_W_ANCOM.csv", row.names = T)
write.csv(resfam$beta, "res_beta_ANCOM.csv", row.names = T)

#overall ANCOMBC results at family level
write.csv(res.global.fam, "res_global_ANCOM.csv", row.names = T)


#### Figure 4 ####
#Import combined ANCOM by life stage
wide <- read.csv("res_W_ANCOM_LS.csv")

# subset where W-stat is 100+
wide_sub <- subset(wide, W_overall > 100)

wide_sub$lcf_groups =  ifelse(wide_sub$lfc_p == TRUE & wide_sub$lfc_m == FALSE & wide_sub$lfc_a == FALSE, "p",
                              ifelse(wide_sub$lfc_p == TRUE & wide_sub$lfc_m == TRUE & wide_sub$lfc_a == FALSE, "p + m",
                                     ifelse(wide_sub$lfc_p == TRUE & wide_sub$lfc_m == TRUE & wide_sub$lfc_a == TRUE, "p + m + a",
                                            ifelse(wide_sub$lfc_p == FALSE & wide_sub$lfc_m == TRUE & wide_sub$lfc_a == FALSE, "m",
                                                   ifelse(wide_sub$lfc_p == FALSE & wide_sub$lfc_m == TRUE & wide_sub$lfc_a == TRUE, "m + a",
                                                          ifelse(wide_sub$lfc_p == FALSE & wide_sub$lfc_m == FALSE & wide_sub$lfc_a == TRUE, "a",
                                                                 ifelse(wide_sub$lfc_p == TRUE & wide_sub$lfc_m == FALSE & wide_sub$lfc_a == TRUE, "p + a",NA)
                                                          ))))))


fig4a = ggplot(data = wide_sub, aes(x = W_overall, y = reorder(taxa,W_overall), color = lcf_groups))+
  geom_segment(aes(x=0, xend=W_overall, y=reorder(taxa,W_overall), yend=reorder(taxa,W_overall))) +
  geom_point()+
  theme_classic()+
  theme(text = element_text(size = 16),
        legend.position = c(0.8,0.2),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))+
  labs( x = "ANCOM W statistic", y = "Bacterial family")+
  guides(color=guide_legend(title="Life stage \ndifferences"))
fig4a

beta_long <- data_long <- melt(wide_sub,
                               # ID variables - all the variables to keep but not split apart on
                               id.vars=c("taxa", "W_overall"),
                               # The source columns
                               measure.vars=c("p_beta","m_beta","a_beta" ),
                               # Name of the destination column that will identify the original
                               # column that the measurement came from
                               variable.name="group",
                               value.name="beta")
beta_long
beta_long$groups = ifelse(beta_long$group == "p_beta","Paedomorph",
                          ifelse(beta_long$group == "m_beta","Metamorphosing",
                                 ifelse(beta_long$group == "a_beta","Adult",NA)))
beta_long$taxgroup = paste(beta_long$taxa, beta_long$groups)

qval_long <- data_long <- melt(wide_sub,
                               # ID variables - all the variables to keep but not split apart on
                               id.vars=c("taxa"),
                               # The source columns
                               measure.vars=c("p_qval","m_qval","a_qval"),
                               # Name of the destination column that will identify the original
                               # column that the measurement came from
                               variable.name="group",
                               value.name="qval")
qval_long
qval_long$groups = ifelse(qval_long$group == "p_qval","Paedomorph",
                          ifelse(qval_long$group == "m_qval","Metamorphosing",
                                 ifelse(qval_long$group == "a_qval","Adult",NA)))
qval_long$taxgroup = paste(qval_long$taxa, qval_long$groups)


long = merge(beta_long[c(1,2,4,5,6)], qval_long[c(3,5)], by = "taxgroup")
long$sig = ifelse(long$qval < 0.05, "*", "")

names(long)
color <- ifelse(long$beta < 0, "#FFCC33", "#003366")

fig4b = ggplot(long, aes(x = reorder(taxa, W_overall), y = beta)) +
  geom_col(aes(fill = beta > 0),
           position = "dodge",
           col = "transparent")+
  geom_text(aes(label = sig,
                hjust = ifelse(beta < 0, 1.25, -0.25),
                vjust = 0.7),
            size = 6) +
  coord_flip() +
  labs(x="", y = "Log fold change") + 
  facet_grid(~groups)+
  scale_y_continuous(limits = c(min(long$beta) - 0.2,
                                max(long$beta) + 0.2))+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_rect(colour=NA, fill=NA),
        text = element_text(size=16),
        axis.text.y = element_blank())+
  scale_fill_manual(values=c("#003366","#FFCC33"))
fig4b

fig4 <- (fig4a + fig4b) + 
  plot_layout(widths = c(1,3))+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=12))#add figure labels
fig4 #view multi-panel figure 

ggsave("Fig4.tiff", units="in", width=14, height=6, dpi=300, compression = 'lzw')

# C'est fini!