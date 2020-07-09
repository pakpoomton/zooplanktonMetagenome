# specify working directory
setwd("/Users/pakpoomton/qiime/zooplankton/otus/")

# name .biom and .tre files
biomName = "otu_table_mc2_w_tax_no_pynast_failures.biom"
treeName = "rep_set.tre"
sampleDataName = "sampleDataZooplankton.csv"

# specify which taxonomic level to create stack plot (Phylum, Class, Order, Family, Genus, Species)
OTUlevel = "Rank3"

# do not show data on stack bar if count (%) of a taxa is blow this CUTOFF
CUTOFF = 0.01

# scaling plot margins
sc=2

relOTU = TRUE

###########################################
####### load libraries ###################

library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(matrixStats)
library(ape)
theme_set(theme_bw())

########################
###### import data #####
biomOTUTax = import_biom(biomName)  # otu table & taxonomic table
phyloTree = read.tree(treeName)     # phylogenetic tree

sampleData = read.csv(sampleDataName, sep = ',', header = T)

sDat = sample_data(sampleData)
row.names(sDat) = sampleData$SampleID
  
# create phyloseq object
physeq = phyloseq(otu_table(biomOTUTax), tax_table(biomOTUTax), phyloTree, sDat)

# group data by pre-specified OTUlevel
phyloGlom = tax_glom(physeq, OTUlevel)

###

## for taxonomic matrix, exclude data at levels below OTUlevel
TAXmt = tax_table(phyloGlom)
indTaxLevel = match(OTUlevel,colnames(TAXmt))
TAXm = TAXmt[, 1:indTaxLevel]
## get otu table (after grouping)
OTUm = otu_table(phyloGlom, taxa_are_rows = TRUE)

# convert absolute to relative count as specified
if (relOTU) {
  OTUm = OTUm / rep(colSums(OTUm), each = nrow(OTUm))
  # screen out taxa with abundance lower than CUTOFF
  idxSel  = rowMaxs(as.matrix(OTUm)) > CUTOFF
  OTUm = OTUm[idxSel, ]
  TAXm = TAXm[idxSel, ]
}

physeqm = phyloseq(OTUm, TAXm) ## phyloseq object after grouping 

# order of bar plot
desired_order = c("G1", "G2", "G3", "R1", "R2", "R3", 
              "S1", "S2", "S3")

######################################################
########### Processing data before plotting #########

set.seed(123458)
colorset = sample(wes_palette(length(rownames(OTUm)), name = "Darjeeling1", type = "continuous"))


# generate abundance bar plots
p <- plot_bar(physeqm, fill = OTUlevel)  + theme(plot.margin = margin(6*sc,1*sc,6*sc,1*sc,"cm")) +
  theme(legend.position="right") + guides(fill=guide_legend(ncol=1)) +
  theme(text = element_text(size=20)) + scale_fill_manual(values = colorset)

# take care of plotting order
pd <- p$data
pd$Sample <- factor(pd$Sample, levels = desired_order)
p$data <- pd
print(p)

ggsave(filename = "myplot.png", plot = last_plot(), 
       width=40, height=40, unit="cm")

###########################

#Remove OTUs that do not show appear more than 2 times in more than 1/3 the samples
wh0 = genefilter_sample(physeq, filterfun_sample(function(x) x > 5), A=0.33*nsamples(physeq))
physeq1 = prune_taxa(wh0, physeq)

# Transform to even sampling depth.
physeq1 = transform_sample_counts(physeq1, function(x) 1E6 * x/sum(x))

# Keep only the most abundant five phyla.
phylum.sum = tapply(taxa_sums(physeq1), tax_table(physeq1)[, "Rank2"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
physeq1 = prune_taxa((tax_table(physeq1)[, "Rank2"] %in% top5phyla), physeq1)

## MDS (“PCoA”) on Unifrac Distances
ordu = ordinate(physeq1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(physeq1, ordu, , color="SampleType")

### ref. https://github.com/joey711/phyloseq/issues/616 
