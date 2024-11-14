###### Script of 16S rRNA data analysis of pig gut microbiota analysis after phyloseq object generation using DADA2 pipeline available at https://benjjneb.github.io/dada2/ ######
##### In this study the data was rarified to 39K reads/sample #####

# Load packages
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2"); theme_set(theme_bw())
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ape)
library(vegan)
library(microbiome)
library(dplyr)
library(Biostrings)
library(viridis)
library(RColorBrewer)
library(knitr)
library(forcats)
library(microbiomeutilities)
library(microbiomeMarker)
library(speedyseq)
library(ggsci)
library(ggplot2)

# Set working directory
setwd("/PATH")

#Load phyloseq object
ps <- readRDS("ps.rds")

### Subset phyloseq
# Based in Timepoint
ps.day_neg7 <- subset_samples(physeq = ps, Timepoint == 'day_-7')
ps.day_0 <- subset_samples(physeq = ps, Timepoint == 'day_0')
ps.day_14 <- subset_samples(physeq = ps, Timepoint == 'day_14')

# Based on Treatment
ps.neg <- subset_samples(physeq = ps, Treatment == 'Neg control')
ps.glu1 <- subset_samples(physeq = ps, Treatment == '1% Glu')
ps.glu2 <- subset_samples(physeq = ps, Treatment == '2% Glu')
ps.asp1 <- subset_samples(physeq = ps, Treatment == '1% Asp')
ps.asp2 <- subset_samples(physeq = ps, Treatment == '2% Asp')
ps.ab <- subset_samples(physeq = ps, Treatment == '1% Antibiotics')
ps.pos <- subset_samples(physeq = ps, Treatment == 'Pos control')

## Alpha diversity
# Define sample comparisons
sampleCompare <- list(c('Neg control','2% Asp'), c('Neg control','1% Glu'), c('Neg control','1% Antibiotics'), c('Neg control','1% Asp'), c('Neg control','2% Glu'), c('Neg control','Pos control'))
sampleCompare2 <- list(c('day_-7','day_0'), c('day_0','day_14'), c('day_-7','day_14'))

# Alpha diversity based on timepoint

plot_richness(ps, x="Timepoint", measures=c("Shannon", "Simpson","Observed", fill="Timepoint"))+geom_boxplot(lwd=1, aes(fill=Timepoint))+ scale_fill_manual(values=c("#009E73", "#0072B2", "#E69F00"))+ stat_compare_means(comparisons = sampleCompare2, method = "wilcox.test", size=4, label="p.signif")+ geom_point(size=1.5, alpha=2.5)+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold") ,strip.text = element_text(size = 12, colour = "black", face = "bold"), legend.title = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 12, face = "bold", colour = "black"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"))



# Alpha diversity based on treatment group at each sampling timepoint

# day -7
p = plot_richness(ps.day_neg7, x="Treatment", measures=c("Shannon", "Simpson","Observed", fill="Treatment"))+geom_boxplot(lwd=1, aes(fill=Treatment))+ scale_fill_manual(values=c("#a6cee3", "#b2df8a", "#fdbf6f","#006600", "#9900CC", "#FF9900", "red"))+ stat_compare_means(comparisons = sampleCompare, method = "wilcox.test", size=4, label="p.signif")+ geom_point(size=1.5, alpha=2.5)+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.2, hjust = 1.0, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold") ,strip.text = element_text(size = 12, colour = "black", face = "bold"), legend.title = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 12, face = "bold", colour = "black"), axis.text.y = element_text(size = 12, colour = "black", face = "bold")) + ggtitle("Alpha diversity Day -7")

newSTorder = c("Neg control", "Pos control", "1% Antibiotics", 
               "1% Glu", "2% Glu", "1% Asp", "2% Asp")

p$data$Treatment <- as.character(p$data$Treatment)
p$data$Treatment <- factor(p$data$Treatment, levels=newSTorder)
print(p)





# day 0
q= plot_richness(ps.day_0, x="Treatment", measures=c("Shannon", "Simpson","Observed", fill="Treatment"))+geom_boxplot(lwd=1, aes(fill=Treatment))+ scale_fill_manual(values=c("#a6cee3", "#b2df8a", "#fdbf6f","#006600", "#9900CC", "#FF9900", "red"))+ stat_compare_means(comparisons = sampleCompare, method = "wilcox.test", size=4, label="p.signif")+ geom_point(size=1.5, alpha=2.5)+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.2, hjust = 1.0, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold") ,strip.text = element_text(size = 12, colour = "black", face = "bold"), legend.title = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 12, face = "bold", colour = "black"), axis.text.y = element_text(size = 12, colour = "black", face = "bold")) + ggtitle("Alpha diversity Day 0")


newSTorder = c("Neg control", "Pos control", "1% Antibiotics", 
               "1% Glu", "2% Glu", "1% Asp", "2% Asp")

q$data$Treatment <- as.character(q$data$Treatment)
q$data$Treatment <- factor(q$data$Treatment, levels=newSTorder)
print(q)





# day 14
r= plot_richness(ps.day_14, x="Treatment", measures=c("Shannon", "Simpson","Observed", fill="Treatment"))+geom_boxplot(lwd=1, aes(fill=Treatment))+ scale_fill_manual(values=c("#a6cee3", "#b2df8a", "#fdbf6f","#006600", "#9900CC", "#FF9900", "red"))+ stat_compare_means(comparisons = sampleCompare, method = "wilcox.test", size=4, label="p.signif")+ geom_point(size=1.5, alpha=2.5)+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.2, hjust = 1.0, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold") ,strip.text = element_text(size = 12, colour = "black", face = "bold"), legend.title = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 12, face = "bold", colour = "black"), axis.text.y = element_text(size = 12, colour = "black", face = "bold")) + ggtitle("Alpha diversity Day 14")


newSTorder = c("Neg control", "Pos control", "1% Antibiotics", 
               "1% Glu", "2% Glu", "1% Asp", "2% Asp")

r$data$Treatment <- as.character(r$data$Treatment)
r$data$Treatment <- factor(r$data$Treatment, levels=newSTorder)
print(r)


## Beta diversity
# PCoA plot (weighted unifrac)
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
plot(random_tree)
physeq1 = merge_phyloseq(ps, random_tree)
physeq1

ordu = ordinate(physeq1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(physeq1, ordu, color="Timepoint") + geom_point(size=2, alpha=0.9) + theme_bw() + stat_ellipse(geom = "polygon", linetype = 2,alpha=0, aes(fill=Timepoint))+ ggtitle("PCoA Plot - weighted unifrac") +  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00")) +  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00")) + theme(axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.2, hjust = 1.0, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_text(size = 12, face = "bold") ,strip.text = element_text(size = 12, colour = "black", face = "bold"), legend.title = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 12, face = "bold", colour = "black"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"))

# Check for significance using Adonis2
samdf <- as(sample_data(ps), "data.frame")

adonis_pp <- adonis2(phyloseq::distance(ps, method="bray") ~ Timepoint,
                     data = samdf)
adonis_pp

# Distance from centroid
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

# export data from phyloseq to vegan-compatible object
pp_vegan <- veganotu(ps)

# make a data frame that can be used in vegan from the sample_data in the
# phyloseq object
sampledf <- data.frame(sample_data(ps))
pp_BC <- vegdist(wisconsin(sqrt(pp_vegan)), method = "bray")
betadisp_pp <- betadisper(pp_BC, sampledf$Timepoint)
betadisp_pp

par(mar = c(7, 4, 2, 2))
boxplot(betadisp_pp, xlab = "", las = 2, cex.axis = 0.8)


#Anova
anova(betadisp_pp)

#Tukey multiple comparisons
TukeyHSD(betadisp_pp)

## Genus level bar plot

# Define functions
setOrder <- function(df,taxaLevel){
  temp <- data.frame(Abundance=df$Abundance,taxaLevel=df[,taxaLevel])
  temp <- temp %>% group_by(taxaLevel) %>% summarise_all(sum)
  temp <- as.data.frame(temp)
  taxaorder <- temp[,1][order(temp[,2],decreasing = T)]
  print(taxaorder)
}

addOther <- function(psObj){
  temp.otu = as.data.frame(t(otu_table(psObj)))
  temp.otu['Other',] = 1 - sample_sums(psObj)
  temp.taxa = as.data.frame(tax_table(psObj))
  temp.taxa['Other',] <- rep('Other',length(rank_names(psObj)))
  temp.ps <- phyloseq(
    otu_table(temp.otu, taxa_are_rows = T),
    tax_table(as.matrix(temp.taxa)),
    sample_data(psObj)
  )
}

# Define colours
colours = c("coral3", "violet","salmon3", "cadetblue3", "burlywood3","brown",
          "skyblue", "cadetblue","goldenrod4", "gold", "chartreuse", "salmon",
          "seagreen", "orange","azure", "cornsilk4", "limegreen", "steelblue", "maroon3", "tomato", "#009E73")

# Get top 20 Genera
temp.taxa <- as.data.frame(tax_table(ps))
temp.taxa$Genus <- gsub(pattern = '_.*',replacement = '',x = temp.taxa$Genus,perl = T)
tax_table(ps) <- tax_table(as.matrix(temp.taxa))
ps.genus <- tax_glom(physeq = ps, taxrank = "Genus")
ps.genus.na <- subset_taxa(ps.genus, Genus!="NA")
ps.genus.rel <- transform_sample_counts(ps.genus.na, function(OTU) OTU/sum(OTU))
top20 <- names(sort(taxa_sums(ps.genus.rel), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps.genus.rel)
# add other
ps.top20.w.other <- addOther(ps.top20)
nb.cols <- length(taxa_names(ps.top20.w.other))
nb.cols
my_color = colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ps.genus.melt <- psmelt(ps.top20.w.other)
ps.genus.filt <- ps.genus.melt[,"Abundance" > 0]

# Plot
p= ggplot(ps.genus.filt, aes(x = PenID, y = Abundance, fill = fct_reorder(Genus, Abundance))) + geom_bar(stat = "identity", colour= "black")+ facet_grid(Timepoint~factor(Treatment, levels=c("Neg control", "Pos control", "1% Antibiotics","1% Glu", "2% Glu", "1% Asp", "2% Asp")), scales="free", space="free_x")+ theme(axis.text.x = element_text(angle = 90, size = 16, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 22, face = "bold"),strip.text = element_text(size = 22, colour = "black", face = "bold"), legend.title = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 18, face = "bold", colour = "black"), axis.text.y = element_text(size = 20, colour = "black", face = "bold"))+ scale_y_continuous()+ labs(x="", y="Relative_Abundance", fill="") +
  theme(legend.position = "bottom") +theme(legend.direction = "horizontal") + scale_fill_manual(values = colours)


newSTorder = c("Neg control", "Pos control", "1% Antibiotics",
               "1% Glu", "2% Glu", "1% Asp", "2% Asp")

p$data$Treatment <- as.character(p$data$Treatment)
p$data$Treatment <- factor(p$data$Treatment, levels=newSTorder)
print(p)


## Differentially abundant taxa

# subset phyloseq based on timepoint
ps.day_14 <- subset_samples(physeq = ps, Timepoint == 'day_14')

ps.day_14@tax_table <- gsub(" ", "_", ps.day_14@tax_table)
ps.day_14@tax_table <- gsub("[", "",ps.day_14@tax_table, fixed = TRUE)
ps.day_14@tax_table <- gsub("]", "",ps.day_14@tax_table, fixed = TRUE)

ps.genus14 <- tax_glom(physeq = ps.day_14, taxrank = "Genus")

tax.tab14<- as.data.frame(ps.genus14@tax_table) %>% select(6)
tax_table(ps.genus14) <- as.matrix(tax.tab5)

# Run Anova
mm_anova_genus <- run_test_multiple_groups(
  ps.genus14,
  group = "Treatment",
  method = "anova"
)

anova_day14_genus <- as.data.frame(mm_anova_genus@marker_table)
write.csv(anova_day14_genus, "anova_day14_genus.csv")

# Run Tukey's pairwise comparison
pht14 <- run_posthoc_test(ps.genus14, group = "Treatment", norm = "TSS", norm_para = list(), conf_level = 0.95, method= c("tukey"))

posthoc14_genus <- as.data.frame(pht14@result)
write.csv(posthoc14_genus, "posthoc_genus_14.csv")

# List of significantly different genera
gen_to_keep_14=c("Alistipes", "Clostridium sensu stricto 1","Candidatus Stoquefichus", "Eisenbergiella","Hungatella","Limosilactobacillus", "Megasphaera","Mitsuokella", "Monoglobus","Rikenellaceae RC9 gut group","Ruminococcus","Streptococcus","Terrisporobacter","UCG-003")

# convert phyloseq to relative abundance
ps.rel <- ps %>% transform_sample_counts(function(x) {x/sum(x)*100})

# Glom at genus level
ps.rel.genus <- tax_glom(physeq = ps.rel, taxrank = "Genus")

# Subset for taxa of interest
rel.g3<- subset_taxa(ps.rel.genus, Genus %in% gen_to_keep_14)

# Plot
y= phyloseq::psmelt(rel.g3) %>% dplyr::select(2:19) %>% relocate(Genus) %>% 
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") + scale_colour_manual(values=c("#a6cee3", "#b2df8a", "#fdbf6f","#006600", "#9900CC", "#FF9900", "red")) +
  facet_wrap(~ Genus, scales = "free_y")+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.2, hjust = 1.0, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), axis.title.x = element_text(size = 12, face = "bold") ,strip.text = element_text(size = 10, colour = "black", face = "bold"), legend.title = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 12, face = "bold", colour = "black"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"))

newSTorder = c("Neg control", "Pos control", "1% Antibiotics", 
               "1% Glu", "2% Glu", "1% Asp", "2% Asp")

y$data$Treatment <- as.character(y$data$Treatment)
y$data$Treatment <- factor(y$data$Treatment, levels=newSTorder)
print(y)

# END
