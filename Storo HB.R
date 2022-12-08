#https://astrobiomike.github.io/amplicon/dada2_workflow_ex#processing-with-dada2-in-r
#set up working environment
library(dada2)
Storo.path <- "/Volumes/BONO/Internship_2022/SRR_DATA/Storo_fastq"
list.files(Storo.path)
Storo.fnFs <- sort(list.files(Storo.path, pattern="_pass.fastq", full.names = TRUE))
Storo.sample.names <- sapply(strsplit(basename(Storo.fnFs), "_"), `[`, 1)
Storo.sample.names
#quality trimming/filtering
plotQualityProfile(Storo.fnFs[1:2])
Storo.filtFs <- file.path(Storo.path, "filtered", paste0(Storo.sample.names, "_filt.fastq.gz"))
Storo.filtered.out <- filterAndTrim(Storo.fnFs, Storo.filtFs,
                              maxEE=c(2),
                              rm.phix=TRUE, minLen=175, truncLen=c(245))
class(Storo.filtered.out)
dim(Storo.filtered.out)
Storo.filtered.out
plotQualityProfile(Storo.filtFs[1:2])
#error model
Storo.errFs <- learnErrors(Storo.filtFs)
plotErrors(Storo.errFs, nominalQ=TRUE)
#dereplication
Storo.derepFs <- derepFastq(Storo.filtFs, verbose=TRUE)
names(Storo.derepFs) <- Storo.sample.names
#inferring ASVs
Storo.dadaFs <- dada(Storo.derepFs, err=Storo.errFs, pool="pseudo")
#count table
Storo.seqtab <- makeSequenceTable(Storo.dadaFs)
class(Storo.seqtab)
dim(Storo.seqtab)
Storo.seqtab
#chimera identification
Storo.seqtab.nochim <- removeBimeraDenovo(Storo.seqtab, verbose=TRUE)
sum(Storo.seqtab.nochim)/sum(Storo.seqtab)
#overview of counts
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
summary_tab
write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)
#assigning taxonomy
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("DECIPHER")
library(DECIPHER)
Storo.taxa <- assignTaxonomy(Storo.seqtab.nochim, "/Volumes/BONO/Internship_2022/SRR_DATA/Jose_fastq/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
Storo.taxa <- addSpecies(Storo.taxa, "/Volumes/BONO/Internship_2022/SRR_DATA/Jose_fastq/silva_species_assignment_v132.fa.gz")
Storo.taxa.print <- Storo.taxa
rownames(Storo.taxa.print) <- NULL
head(Storo.taxa.print)
Storo.taxa.print
#extracting goods from DADA2 to other formats
Storo.asv.seqs <- colnames(Storo.seqtab.nochim)
Storo.asv.headers <- vector(dim(Storo.seqtab.nochim)[2], mode="character")
for (i in 1:dim(Storo.seqtab.nochim)[2]) {Storo.asv.headers[i] <- paste(">ASV", i, sep="_")}
#fasta file extract
Storo.asv.fasta <- c(rbind(Storo.asv.headers, Storo.asv.seqs))
write(Storo.asv.fasta, "Storo.ASVs.fa")
#count table extract
Storo.asv.tab <- t(Storo.seqtab.nochim)
row.names(Storo.asv.tab) <- sub(">", "", Storo.asv.headers)
write.table(Storo.asv.tab, "Storo.ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
#tax table and setting unclassified
Storo.ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

write.matrix(Storo.taxa, file="ASVs_taxonomy.tsv")
write.matrix(Storo.taxa, file="ASVs_Taxa(raw).csv")
installed.packages("MASS")
#re reading in data
library(tidyverse) ; packageVersion("tidyverse")
library(phyloseq) ; packageVersion("phyloseq")
library(vegan) ; packageVersion("vegan")
BiocManager::install("DESeq2")
library(DESeq2) ; packageVersion("DESeq2")
library(dendextend) ; packageVersion("dendextend")
library(viridis) ; packageVersion("viridis")

Storo.count.tab <- read.table("Storo.ASVs_counts.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")
Storo.tax.tab <- as.matrix(read.csv("Storo.ASVs_Taxa.csv", header=T,
                                      check.names=F))
#sample_info_tab is created from study metadata
Storo.sample.info.tab <- read.table("Storo_Sample_Info.csv", header=T, row.names=1,
                              check.names=F, sep="\t")
Storo.sample.info.tab$Color <- as.character(Storo.sample.info.tab$Color)
dim(Storo.sample.info.tab)
Storo.sample.info.tab

#BETA DIVERSITY
#normalizing for sampling depth
Storo.deseq.counts <- DESeqDataSetFromMatrix(Storo.count.tab, colData = Storo.sample.info.tab, design = ~Host)
Storo.deseq.counts.vst <- varianceStabilizingTransformation(Storo.deseq.counts)
  #fixing error
Storo.deseq.counts <- estimateSizeFactors(Storo.deseq.counts, type = "poscounts")
Storo.deseq.counts.vst <- varianceStabilizingTransformation(Storo.deseq.counts)
Storo.vst.trans.count.tab <- assay(Storo.deseq.counts.vst)
Storo.euc.dist <- dist(t(Storo.vst.trans.count.tab))
#hierarchical clustering
Storo.euc.clust <- hclust(Storo.euc.dist, method="ward.D2")
plot(Storo.euc.clust) 
Storo.euc.dend <- as.dendrogram(Storo.euc.clust, hang=0.1)
Storo.dend.cols <- as.character(Storo.sample.info.tab$Color[order.dendrogram(Storo.euc.dend)])
labels_colors(Storo.euc.dend) <- Storo.dend.cols
plot(Storo.euc.dend, ylab="VST Euc. dist.")
#ordination
Storo.vst.count.phy <- otu_table(Storo.vst.trans.count.tab, taxa_are_rows=T)
Storo.sample.info.tab.phy <- sample_data(Storo.sample.info.tab)
Storo.vst.physeq <- phyloseq(Storo.vst.count.phy, Storo.sample.info.tab.phy)
Storo.vst.pcoa <- ordinate(Storo.vst.physeq, method="MDS", distance="euclidean")
Storo.eigen.vals <- Storo.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.vst.physeq, Storo.vst.pcoa, color="Host") + 
  geom_jitter(aes(color=factor(Host))) + geom_point(size=1) + labs(col="Host") + 
  coord_fixed(sqrt(Storo.eigen.vals[2]/Storo.eigen.vals[1])) + ggtitle("PCoA by Species and Environment") + 
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "top")

#ALPHA DIVERSITY
rarecurve(t(Storo.count.tab), step=100, col=Storo.sample.info.tab$Color, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(t(Storo.count.tab)))))
#richness and diversity estimates
Storo.count.tab.phy <- otu_table(Storo.count.tab, taxa_are_rows=T)
Storo.tax.tab.phy <- tax_table(Storo.tax.tab) #error correction next
rownames(Storo.tax.tab) <- paste("ASV", 1:nrow(Storo.count.tab), sep="_")
Storo.ASV.physeq <- phyloseq(Storo.count.tab.phy, Storo.tax.tab.phy, Storo.sample.info.tab.phy)
plot_richness(Storo.ASV.physeq, x="Host", measures=c("Chao1", "Shannon")) + geom_boxplot() +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plot_richness(Storo.ASV.physeq, x="Isolation_Source", color="Host", measures=c("Chao1", "Shannon")) +
  scale_color_manual(values=unique(Storo.sample.info.tab$Color[order(Storo.sample.info.tab$Host)])) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

anova(betadisper(Storo.euc.dist, Storo.sample.info.tab$Organism))

#taxonomic summaries (phyla based on example, see replicates file for genus and species)
Storo.phyla.counts.tab <- otu_table(tax_glom(Storo.ASV.physeq, taxrank="Phylum")) 
Storo.phyla.tax.vec <- as.vector(tax_table(tax_glom(Storo.ASV.physeq, taxrank="Phylum"))[,"Phylum"]) 
rownames(Storo.phyla.counts.tab) <- as.vector(Storo.phyla.tax.vec)
Storo.unclassified.tax.counts <- colSums(Storo.count.tab) - colSums(Storo.phyla.counts.tab)
Storo.phyla.and.unidentified.counts.tab <- rbind(Storo.phyla.counts.tab, "Unclassified"=Storo.unclassified.tax.counts)
Storo.temp.major.taxa.counts.tab <- Storo.phyla.and.unidentified.counts.tab[!row.names(Storo.phyla.and.unidentified.counts.tab) %in% "Proteobacteria", ]
Storo.class.counts.tab <- otu_table(tax_glom(Storo.ASV.physeq, taxrank="Class"))
Storo.class.tax.phy.tab <- tax_table(tax_glom(Storo.ASV.physeq, taxrank="Class"))
Storo.phy.tmp.vec <- Storo.class.tax.phy.tab[,2]
Storo.class.tmp.vec <- Storo.class.tax.phy.tab[,3]
Storo.rows.tmp <- row.names(Storo.class.tax.phy.tab)
Storo.class.tax.tab <- data.frame("Phylum"=Storo.phy.tmp.vec, "Class"=Storo.class.tmp.vec, row.names = Storo.rows.tmp)
Storo.proteo.classes.vec <- as.vector(Storo.class.tax.tab[Storo.class.tax.tab$Phylum == "Proteobacteria", "Class"])
rownames(Storo.class.counts.tab) <- as.vector(Storo.class.tax.tab$Class) 
Storo.proteo.class.counts.tab <- Storo.class.counts.tab[row.names(Storo.class.counts.tab) %in% Storo.proteo.classes.vec, ]
Storo.proteo.no.class.annotated.counts <- Storo.phyla.and.unidentified.counts.tab[row.names(Storo.phyla.and.unidentified.counts.tab) %in% "Proteobacteria", ] - colSums(Storo.proteo.class.counts.tab)
Storo.major.taxa.counts.tab <- rbind(Storo.temp.major.taxa.counts.tab, Storo.proteo.class.counts.tab, "Unresolved_Proteobacteria"=Storo.proteo.no.class.annotated.counts)
identical(colSums(Storo.major.taxa.counts.tab), colSums(Storo.count.tab)) 
Storo.major.taxa.proportions.tab <- apply(Storo.major.taxa.counts.tab, 2, function(x) x/sum(x)*100)
dim(Storo.major.taxa.proportions.tab)
Storo.temp.filt.major.taxa.proportions.tab <- data.frame(Storo.major.taxa.proportions.tab[apply(Storo.major.taxa.proportions.tab, 1, max) > 5, ])
dim(Storo.temp.filt.major.taxa.proportions.tab) 
Storo.filtered.proportions <- colSums(Storo.major.taxa.proportions.tab) - colSums(Storo.temp.filt.major.taxa.proportions.tab)
Storo.filt.major.taxa.proportions.tab <- rbind(Storo.temp.filt.major.taxa.proportions.tab, "Other"=Storo.filtered.proportions)
Storo.filt.major.taxa.proportions.tab.for.plot <- Storo.filt.major.taxa.proportions.tab
Storo.filt.major.taxa.proportions.tab.for.plot$Major_Taxa <- row.names(Storo.filt.major.taxa.proportions.tab.for.plot)
Storo.filt.major.taxa.proportions.tab.for.plot.g <- pivot_longer(Storo.filt.major.taxa.proportions.tab.for.plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()
head(Storo.filt.major.taxa.proportions.tab.for.plot.g)
head(Storo.filt.major.taxa.proportions.tab.for.plot)
Storo.sample.info.for.merge<-data.frame("Sample"=row.names(Storo.sample.info.tab), "Host"=Storo.sample.info.tab$Host, "Color"=Storo.sample.info.tab$Color, "Isolation_Source"=Storo.sample.info.tab$Isolation_Source, stringsAsFactors=F)
Storo.filt.major.taxa.proportions.tab.for.plot.g2 <- merge(Storo.filt.major.taxa.proportions.tab.for.plot.g, Storo.sample.info.for.merge)
ggplot(Storo.filt.major.taxa.proportions.tab.for.plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="All samples")
ggplot(Storo.filt.major.taxa.proportions.tab.for.plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(Host), shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.tab.for.plot.g2$Color[order(Storo.filt.major.taxa.proportions.tab.for.plot.g2$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", title="Phylum Taxonomy Across All Samples")

#SPECIES ONLY PLOTS#
#Negaprion brevirostris samples
Storo.NB.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Negaprion_brevirostris"]
Storo.filt.major.taxa.proportions.NB.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.NB.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.NB.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.NB.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.NB.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Negaprion brevirostris samples only")
ggplot(Storo.filt.major.taxa.proportions.NB.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Negaprion brevirostris samples only")
#Carcharhinus melanopterus samples
Storo.CM.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Carcharhinus_melanopterus"]
Storo.filt.major.taxa.proportions.CM.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CM.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.CM.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.CM.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.CM.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Carcharhinus melanopterus samples only")
ggplot(Storo.filt.major.taxa.proportions.CM.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus melanopterus samples only")
#Carcharhinus leucas samples
Storo.CL.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Carcharhinus_leucas"]
Storo.filt.major.taxa.proportions.CL.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CL.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.CL.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.CL.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.CL.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Carcharhinus leucas samples only")
ggplot(Storo.filt.major.taxa.proportions.CL.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus leucas samples only")
#Carcharhinus perezii samples
Storo.CPe.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Carcharhinus_perezii"]
Storo.filt.major.taxa.proportions.CPe.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CPe.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.CPe.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.CPe.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.CPe.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Carcharhinus perezii samples only")
ggplot(Storo.filt.major.taxa.proportions.CPe.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus perezii samples only")
#Carcharhinus plumbeus samples
Storo.CPl.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Carcharhinus_plumbeus"]
Storo.filt.major.taxa.proportions.CPl.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CPl.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.CPl.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.CPl.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.CPl.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Carcharhinus plumbeuss samples only")
ggplot(Storo.filt.major.taxa.proportions.CPl.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus plumbeus samples only")
#Galeocerdo cuvier samples, error
Storo.GCu.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Galeocerdo_cuvier"]
Storo.filt.major.taxa.proportions.GCu.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.GCu.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.GCu.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.GCu.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.GCu.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Galeocerdo cuvier samples only")
ggplot(Storo.filt.major.taxa.proportions.GCu.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Galeocerdo cuvier samples only")
#Ginglymostoma cirratum samples
Storo.GCi.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Ginglymostoma_cirratum"]
Storo.filt.major.taxa.proportions.GCi.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.GCi.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.GCi.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.GCi.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.GCi.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Ginglymostoma cirratum samples only")
ggplot(Storo.filt.major.taxa.proportions.GCi.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Ginglymostoma cirratum samples only")
#Rhizoprionodon terraenovae samples
Storo.RT.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Rhizoprionodon_terraenovae"]
Storo.filt.major.taxa.proportions.RT.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.RT.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.RT.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.RT.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.RT.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Rhizoprionodon terraenovae samples only")
ggplot(Storo.filt.major.taxa.proportions.RT.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Rhizoprionodon terraenovae samples only")
#Sarda sarda samples
Storo.SS.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Host == "Sarda_sarda"]
Storo.filt.major.taxa.proportions.SS.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.SS.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.SS.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.SS.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.SS.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Sarda sarda samples only")
ggplot(Storo.filt.major.taxa.proportions.SS.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Sarda sarda samples only")
#LOCATION PLOTS ONLY#
#gills
Storo.gills.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Isolation_Source == "gills"]
Storo.filt.major.taxa.proportions.gills.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.gills.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.gills.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.gills.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.gills.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Gill samples only")
ggplot(Storo.filt.major.taxa.proportions.gills.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Gill samples only")
#cloaca
Storo.cloaca.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Isolation_Source == "cloaca"]
Storo.filt.major.taxa.proportions.cloaca.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.cloaca.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.cloaca.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.cloaca.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.cloaca.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Cloaca samples only")
ggplot(Storo.filt.major.taxa.proportions.cloaca.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Cloaca samples only")
#skin
Storo.skin.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Isolation_Source == "skin"]
Storo.filt.major.taxa.proportions.skin.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.skin.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.skin.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.skin.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.skin.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Skin samples only")
ggplot(Storo.filt.major.taxa.proportions.skin.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Skin samples only")
#teeth
Storo.teeth.sample.IDs <- row.names(Storo.sample.info.tab)[Storo.sample.info.tab$Isolation_Source == "Teeth"]
Storo.filt.major.taxa.proportions.teeth.only.tab.for.plot.g <- Storo.filt.major.taxa.proportions.tab.for.plot.g2[Storo.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.teeth.sample.IDs, ]
ggplot(Storo.filt.major.taxa.proportions.teeth.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.filt.major.taxa.proportions.teeth.only.tab.for.plot.g$Color[order(Storo.filt.major.taxa.proportions.teeth.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Taxa", title="Teeth samples only")
ggplot(Storo.filt.major.taxa.proportions.teeth.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Teeth samples only")

##BETADISPER AND PERMUTATIONAL ANOVA##
#Testing in non-replicates: ignore through remainder of file and refer to replicates file
anova(betadisper(Storo.euc.dist, Storo.sample.info.tab$Organism))

#ISOLATION SOURCE GROUPS#

anova(betadisper(Storo.euc.dist, Storo.sample.info.tab$Isolation_Source))
#gills
Storo.gills.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.gills.sample.IDs]))
Storo.gills.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.gills.sample.IDs, ]
anova(betadisper(Storo.gills.euc.dist, Storo.gills.sample.info.tab$Host))
adonis2(Storo.gills.euc.dist~Storo.gills.sample.info.tab$Host)
Storo.gills.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.gills.sample.IDs], taxa_are_rows=T)
Storo.gills.sample.info.tab.phy <- sample_data(Storo.gills.sample.info.tab)
Storo.gills.vst.physeq <- phyloseq(Storo.gills.vst.count.phy, Storo.gills.sample.info.tab.phy)
Storo.gills.vst.pcoa <- ordinate(Storo.gills.vst.physeq, method="MDS", distance="euclidean")
Storo.gills.eigen.vals <- Storo.gills.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.gills.vst.physeq, Storo.gills.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.gills.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Gills composition between host species") +
  annotate("text", x=300, y=40, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.gills.eigen.vals[2]/Storo.gills.eigen.vals[1])) + ggtitle("PCoA - Gill samples only") + 
  scale_color_manual(values=unique(Storo.gills.sample.info.tab$Color[order(Storo.gills.sample.info.tab$Host)])) + 
  theme_bw() + theme(legend.position="top")
#cloaca
Storo.cloaca.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.cloaca.sample.IDs]))
Storo.cloaca.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.cloaca.sample.IDs, ]
anova(betadisper(Storo.cloaca.euc.dist, Storo.cloaca.sample.info.tab$Host))
adonis2(Storo.cloaca.euc.dist~Storo.cloaca.sample.info.tab$Host) 
Storo.cloaca.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.cloaca.sample.IDs], taxa_are_rows=T)
Storo.cloaca.sample.info.tab.phy <- sample_data(Storo.cloaca.sample.info.tab)
Storo.cloaca.vst.physeq <- phyloseq(Storo.cloaca.vst.count.phy, Storo.cloaca.sample.info.tab.phy)
Storo.cloaca.vst.pcoa <- ordinate(Storo.cloaca.vst.physeq, method="MDS", distance="euclidean")
Storo.cloaca.eigen.vals <- Storo.cloaca.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.cloaca.vst.physeq, Storo.cloaca.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.cloaca.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Cloaca composition between host species") +
  annotate("text", x=300, y=40, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.cloaca.eigen.vals[2]/Storo.cloaca.eigen.vals[1])) + ggtitle("PCoA - Cloaca samples only") + 
  scale_color_manual(values=unique(Storo.cloaca.sample.info.tab$Color[order(Storo.cloaca.sample.info.tab$Host)])) + 
  theme_bw() + theme(legend.position="top")
#skin
Storo.skin.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.skin.sample.IDs]))
Storo.skin.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.skin.sample.IDs, ]
anova(betadisper(Storo.skin.euc.dist, Storo.skin.sample.info.tab$Host))
adonis2(Storo.skin.euc.dist~Storo.skin.sample.info.tab$Host)
Storo.skin.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.skin.sample.IDs], taxa_are_rows=T)
Storo.skin.sample.info.tab.phy <- sample_data(Storo.skin.sample.info.tab)
Storo.skin.vst.physeq <- phyloseq(Storo.skin.vst.count.phy, Storo.skin.sample.info.tab.phy)
Storo.skin.vst.pcoa <- ordinate(Storo.skin.vst.physeq, method="MDS", distance="euclidean")
Storo.skin.eigen.vals <- Storo.skin.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.skin.vst.physeq, Storo.skin.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.skin.eigen.vals[2]/Storo.skin.eigen.vals[1])) + ggtitle("PCoA - Skin samples only") + 
  scale_color_manual(values=unique(Storo.skin.sample.info.tab$Color[order(Storo.skin.sample.info.tab$Host)])) + 
  theme_bw() + theme(legend.position="top")
#teeth
Storo.teeth.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.teeth.sample.IDs]))
Storo.teeth.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.teeth.sample.IDs, ]
anova(betadisper(Storo.teeth.euc.dist, Storo.teeth.sample.info.tab$Host))
adonis2(Storo.teeth.euc.dist~Storo.teeth.sample.info.tab$Host)
Storo.teeth.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.teeth.sample.IDs], taxa_are_rows=T)
Storo.teeth.sample.info.tab.phy <- sample_data(Storo.teeth.sample.info.tab)
Storo.teeth.vst.physeq <- phyloseq(Storo.teeth.vst.count.phy, Storo.teeth.sample.info.tab.phy)
Storo.teeth.vst.pcoa <- ordinate(Storo.teeth.vst.physeq, method="MDS", distance="euclidean")
Storo.teeth.eigen.vals <- Storo.teeth.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.teeth.vst.physeq, Storo.teeth.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.teeth.eigen.vals[2]/Storo.teeth.eigen.vals[1])) + ggtitle("PCoA - Teeth samples only") + 
  scale_color_manual(values=unique(Storo.teeth.sample.info.tab$Color[order(Storo.teeth.sample.info.tab$Host)])) + 
  theme_bw() + theme(legend.position="top")

#HOST GROUPS#
anova(betadisper(Storo.euc.dist, Storo.sample.info.tab$Host))
#Negaprion brevirostris
Storo.NB.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.NB.sample.IDs]))
Storo.NB.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.NB.sample.IDs, ]
anova(betadisper(Storo.NB.euc.dist, Storo.NB.sample.info.tab$Isolation_Source))
adonis(Storo.NB.euc.dist~Storo.NB.sample.info.tab$Isolation_Source) #output does not provide significance value (can't see it)
Storo.NB.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.NB.sample.IDs], taxa_are_rows=T)
Storo.NB.sample.info.tab.phy <- sample_data(Storo.NB.sample.info.tab)
Storo.NB.vst.physeq <- phyloseq(Storo.NB.vst.count.phy, Storo.NB.sample.info.tab.phy)
Storo.NB.vst.pcoa <- ordinate(Storo.NB.vst.physeq, method="MDS", distance="euclidean")
Storo.NB.eigen.vals <- Storo.NB.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.NB.vst.physeq, Storo.NB.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.NB.eigen.vals[2]/Storo.NB.eigen.vals[1])) + ggtitle("PCoA - Negaprion brevirostris samples only") + 
  theme_bw() + theme(legend.position="top")
#Carcharhinus melanopterus samples
Storo.CM.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.CM.sample.IDs]))
Storo.CM.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.CM.sample.IDs, ]
anova(betadisper(Storo.CM.euc.dist, Storo.CM.sample.info.tab$Isolation_Source))
adonis(Storo.CM.euc.dist~Storo.CM.sample.info.tab$Isolation_Source)
Storo.CM.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.CM.sample.IDs], taxa_are_rows=T)
Storo.CM.sample.info.tab.phy <- sample_data(Storo.CM.sample.info.tab)
Storo.CM.vst.physeq <- phyloseq(Storo.CM.vst.count.phy, Storo.CM.sample.info.tab.phy)
Storo.CM.vst.pcoa <- ordinate(Storo.CM.vst.physeq, method="MDS", distance="euclidean")
Storo.CM.eigen.vals <- Storo.CM.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.CM.vst.physeq, Storo.CM.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.CM.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Carcharhinus melanopterus anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.CM.eigen.vals[2]/Storo.CM.eigen.vals[1])) + ggtitle("PCoA - Carcharhinus melanopterus samples only") + 
  theme_bw() + theme(legend.position="top")
#Carcharhinus leucas samples
Storo.CL.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.CL.sample.IDs]))
Storo.CL.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.CL.sample.IDs, ]
anova(betadisper(Storo.CL.euc.dist, Storo.CL.sample.info.tab$Isolation_Source))
adonis(Storo.CL.euc.dist~Storo.CL.sample.info.tab$Isolation_Source) 
Storo.CL.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.CL.sample.IDs], taxa_are_rows=T)
Storo.CL.sample.info.tab.phy <- sample_data(Storo.CL.sample.info.tab)
Storo.CL.vst.physeq <- phyloseq(Storo.CL.vst.count.phy, Storo.CL.sample.info.tab.phy)
Storo.CL.vst.pcoa <- ordinate(Storo.CL.vst.physeq, method="MDS", distance="euclidean")
Storo.CL.eigen.vals <- Storo.CL.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.CL.vst.physeq, Storo.CL.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.CL.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Carcharhinus leucas anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.CL.eigen.vals[2]/Storo.CL.eigen.vals[1])) + ggtitle("PCoA - Carcharhinus leucas samples only") + 
  theme_bw() + theme(legend.position="top")
#Carcharhinus perezii samples
Storo.CPe.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.CPe.sample.IDs]))
Storo.CPe.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.CPe.sample.IDs, ]
anova(betadisper(Storo.CPe.euc.dist, Storo.CPe.sample.info.tab$Isolation_Source))
adonis(Storo.CPe.euc.dist~Storo.CPe.sample.info.tab$Isolation_Source) 
Storo.CPe.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.CPe.sample.IDs], taxa_are_rows=T)
Storo.CPe.sample.info.tab.phy <- sample_data(Storo.CPe.sample.info.tab)
Storo.CPe.vst.physeq <- phyloseq(Storo.CPe.vst.count.phy, Storo.CPe.sample.info.tab.phy)
Storo.CPe.vst.pcoa <- ordinate(Storo.CPe.vst.physeq, method="MDS", distance="euclidean")
Storo.CPe.eigen.vals <- Storo.CPe.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.CPe.vst.physeq, Storo.CPe.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.CPe.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Carcharhinus perezii anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.CPe.eigen.vals[2]/Storo.CPe.eigen.vals[1])) + ggtitle("PCoA - Carcharhinus perezii samples only") + 
  theme_bw() + theme(legend.position="top")
#Carcharhinus plumbeus samples
Storo.CPl.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.CPl.sample.IDs]))
Storo.CPl.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.CPl.sample.IDs, ]
anova(betadisper(Storo.CPl.euc.dist, Storo.CPl.sample.info.tab$Isolation_Source))
adonis(Storo.CPl.euc.dist~Storo.CPl.sample.info.tab$Isolation_Source) 
Storo.CPl.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.CPl.sample.IDs], taxa_are_rows=T)
Storo.CPl.sample.info.tab.phy <- sample_data(Storo.CPl.sample.info.tab)
Storo.CPl.vst.physeq <- phyloseq(Storo.CPl.vst.count.phy, Storo.CPl.sample.info.tab.phy)
Storo.CPl.vst.pcoa <- ordinate(Storo.CPl.vst.physeq, method="MDS", distance="euclidean")
Storo.CPl.eigen.vals <- Storo.CPl.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.CPl.vst.physeq, Storo.CPl.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.CPl.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Carcharhinus plumbeus anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.CPl.eigen.vals[2]/Storo.CPl.eigen.vals[1])) + ggtitle("PCoA - Carcharhinus plumbeus samples only") + 
  theme_bw() + theme(legend.position="top")
#Galeocerdo cuvier samples
Storo.GCu.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.GCu.sample.IDs]))
Storo.GCu.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.GCu.sample.IDs, ]
anova(betadisper(Storo.GCu.euc.dist, Storo.GCu.sample.info.tab$Isolation_Source))
adonis(Storo.GCu.euc.dist~Storo.GCu.sample.info.tab$Isolation_Source) 
Storo.GCu.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.GCu.sample.IDs], taxa_are_rows=T)
Storo.GCu.sample.info.tab.phy <- sample_data(Storo.GCu.sample.info.tab)
Storo.GCu.vst.physeq <- phyloseq(Storo.GCu.vst.count.phy, Storo.GCu.sample.info.tab.phy)
Storo.GCu.vst.pcoa <- ordinate(Storo.GCu.vst.physeq, method="MDS", distance="euclidean")
Storo.GCu.eigen.vals <- Storo.GCu.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.GCu.vst.physeq, Storo.GCu.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.GCu.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Galeocerdo cuvier anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.GCu.eigen.vals[2]/Storo.GCu.eigen.vals[1])) + ggtitle("PCoA - Galeocerdo cuvier samples only") + 
  theme_bw() + theme(legend.position="top")
#Ginglymostoma cirratum samples
Storo.GCi.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.GCi.sample.IDs]))
Storo.GCi.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.GCi.sample.IDs, ]
anova(betadisper(Storo.GCi.euc.dist, Storo.GCi.sample.info.tab$Isolation_Source))
adonis(Storo.GCi.euc.dist~Storo.GCi.sample.info.tab$Isolation_Source)
Storo.GCi.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.GCi.sample.IDs], taxa_are_rows=T)
Storo.GCi.sample.info.tab.phy <- sample_data(Storo.GCi.sample.info.tab)
Storo.GCi.vst.physeq <- phyloseq(Storo.GCi.vst.count.phy, Storo.GCi.sample.info.tab.phy)
Storo.GCi.vst.pcoa <- ordinate(Storo.GCi.vst.physeq, method="MDS", distance="euclidean")
Storo.GCi.eigen.vals <- Storo.GCi.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.GCi.vst.physeq, Storo.GCi.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.GCi.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Ginglymostoma cirratum anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.GCi.eigen.vals[2]/Storo.GCi.eigen.vals[1])) + ggtitle("PCoA - Ginglymostoma cirratum samples only") + 
  theme_bw() + theme(legend.position="top")
#Rhizoprionodon terraenovae samples
Storo.RT.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.RT.sample.IDs]))
Storo.RT.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.RT.sample.IDs, ]
anova(betadisper(Storo.RT.euc.dist, Storo.RT.sample.info.tab$Isolation_Source))
adonis(Storo.RT.euc.dist~Storo.RT.sample.info.tab$Isolation_Source)
Storo.RT.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.RT.sample.IDs], taxa_are_rows=T)
Storo.RT.sample.info.tab.phy <- sample_data(Storo.RT.sample.info.tab)
Storo.RT.vst.physeq <- phyloseq(Storo.RT.vst.count.phy, Storo.RT.sample.info.tab.phy)
Storo.RT.vst.pcoa <- ordinate(Storo.RT.vst.physeq, method="MDS", distance="euclidean")
Storo.RT.eigen.vals <- Storo.RT.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.RT.vst.physeq, Storo.RT.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.RT.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Rhizoprionodon terraenovae anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.RT.eigen.vals[2]/Storo.RT.eigen.vals[1])) + ggtitle("PCoA - Rhizoprionodon terraenovae samples only") + 
  theme_bw() + theme(legend.position="top")
#Sarda sarda samples
Storo.SS.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.SS.sample.IDs]))
Storo.SS.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.SS.sample.IDs, ]
anova(betadisper(Storo.SS.euc.dist, Storo.SS.sample.info.tab$Isolation_Source))
adonis(Storo.SS.euc.dist~Storo.SS.sample.info.tab$Isolation_Source)
Storo.SS.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.SS.sample.IDs], taxa_are_rows=T)
Storo.SS.sample.info.tab.phy <- sample_data(Storo.SS.sample.info.tab)
Storo.SS.vst.physeq <- phyloseq(Storo.SS.vst.count.phy, Storo.SS.sample.info.tab.phy)
Storo.SS.vst.pcoa <- ordinate(Storo.SS.vst.physeq, method="MDS", distance="euclidean")
Storo.SS.eigen.vals <- Storo.SS.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.SS.vst.physeq, Storo.SS.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  geom_text(aes(label=rownames(Storo.SS.sample.info.tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=300, y=68, label="Sarda sarda anatomical location samples") +
  annotate("text", x=300, y=50, label="Permutational ANOVA = ____") + 
  coord_fixed(sqrt(Storo.SS.eigen.vals[2]/Storo.SS.eigen.vals[1])) + ggtitle("PCoA - Sarda sarda samples only") + 
  theme_bw() + theme(legend.position="top")