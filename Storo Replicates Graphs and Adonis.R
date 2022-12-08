#all replicates
Storo.replicate.sample.IDs.tab <- Storo.sample.info.tab[!Storo.sample.info.tab$Host == "Rhizoprionodon_terraenovae" & 
                                                          !Storo.sample.info.tab$Host == "Sarda_sarda" & 
                                                          !Storo.sample.info.tab$Host == "Carcharhinus_melanopterus" & 
                                                          !Storo.sample.info.tab$Host == "Carcharhinus_leucas" &
                                                          !Storo.sample.info.tab$Host == "Seawater",]
Storo.genus.filt.major.taxa.proportions.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.genus.filt.major.taxa.proportions.replicates.only.tab.for.plot.g$Color[order(Storo.genus.filt.major.taxa.proportions.replicates.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(shape=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.species.filt.major.taxa.proportions.replicates.only.tab.for.plot.g$Color[order(Storo.species.filt.major.taxa.proportions.replicates.only.tab.for.plot.g$Host)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Replicate samples only")

#playing with graphs
ggplot(Storo.genus.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(5,100)) + 
  geom_jitter(aes(shape=factor(Host)), size=2, width=0.15, height=0) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Genus Taxa Representation by Shark Host")
ggplot(Storo.genus.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  facet_wrap(~Host, scales = "free_x", nrow = 1) +
  theme(axis.text.x=element_text(angle=90)) + 
  theme(legend.position = "top") +
  theme(legend.text = element_text(size = 8)) +
  theme(legend.key.size = unit(3, "mm")) +
  theme(legend.margin = margin(0, 0, 0, 0, "cm"))+
  labs(x="Sample", title="Major Genus Taxa by Shark Host", fill="Genus")
ggplot(Storo.species.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(5,100)) + 
  geom_jitter(aes(shape=factor(Host)), size=2, width=0.15, height=0) +
  scale_x_discrete(labels=c("atrinae"="Endoziocomonas atrinae", "butzleri"="Acrobacter butzleri", "celer" = "Psychrobacter celer", "clevelandensis" = "Lawsonella clevelandensis", "damselae" = "Photobacterium damselae",
                            "maltophilia" = "Stenotrophomonas maltophilia", "marincola" = "Psychrobacter marincola", "marinus" = "Prochlorococcus MIT9313 marinus", "marisflavi" = "Cohaesibacter marisflavi", "Other" = "Other",
                            "phthalicicus" = "Tropicibacter phthalicicus", "sanguinegens" = "Sneathia sanguinegens", "sanguinis" = "Psychrobacter sanguinis", "sanxanigenens" = "Sphingomonas sanxanigenens","thermosphacta" = "Brochothrix thermosphacta",
                            "Unclassified" = "Unclassified", "vaginalis" = "Gardnerella vaginalis"))+
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="ASV Taxa Representation by Shark Host")
ggplot(Storo.species.filt.major.taxa.proportions.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  facet_wrap(~Host, scales = "free_x", nrow = 1) +
  theme(axis.text.x=element_text(angle=90)) + 
  theme(legend.position = "top", legend.box = "vertical") +
  theme(legend.text = element_text(size = 8)) +
  scale_fill_discrete(labels = c("Endoziocomonas atrinae", "Acrobacter butzleri", "Psychrobacter celer", "Lawsonella clevelandensis","Photobacterium damselae",
                                 "Stenotrophomonas maltophilia","Psychrobacter marincola","Prochlorococcus MIT9313 marinus","Cohaesibacter marisflavi","Other",
                                 "Tropicibacter phtoalicicus","Sneathia sanguinegens","Psychrobacter sanguinis","Sphingomonas sanxanigenens","Brochothrix thermosphacta",
                                 "Unclassified","Gardnerella vaginalis")) +
  labs(x="Sample", title="Major ASV Taxa by Shark Host", fill="ASV Species")


##ISOLATION SOURCES##
anova(betadisper(Storo.euc.dist, Storo.replicate.sample.IDs.tab$Isolation_Source))

#Cloaca
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.cloaca.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.cloaca.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host != "Rhizoprionodon_terraenovae" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Sarda_sarda" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_melanopterus" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_leucas" &
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Seawater" &
                                                                                 Storo.replicate.sample.IDs.tab$Isolation_Source == "cloaca"]
Storo.genus.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.cloaca.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.cloaca.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.genus.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g$Color[order(Storo.genus.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Cloaca replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Cloaca replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.species.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g$Color[order(Storo.species.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Cloaca replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.cloaca.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Cloaca replicate samples only")

Storo.rep.cloaca.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.cloaca.replicate.sample.IDs]))
Storo.rep.cloaca.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.cloaca.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.cloaca.euc.dist, Storo.rep.cloaca.sample.info.tab$Host))
adonis2(Storo.rep.cloaca.euc.dist~Storo.rep.cloaca.sample.info.tab$Host)
Storo.rep.cloaca.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.cloaca.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.cloaca.sample.info.tab.phy <- sample_data(Storo.rep.cloaca.sample.info.tab)
Storo.rep.cloaca.vst.physeq <- phyloseq(Storo.rep.cloaca.vst.count.phy, Storo.rep.cloaca.sample.info.tab.phy)
Storo.rep.cloaca.vst.pcoa <- ordinate(Storo.rep.cloaca.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.cloaca.eigen.vals <- Storo.rep.cloaca.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.cloaca.vst.physeq, Storo.rep.cloaca.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.rep.cloaca.eigen.vals[2]/Storo.rep.cloaca.eigen.vals[1])) + ggtitle("PCoA - Cloaca replicate samples only") + 
  scale_color_manual(values=unique(Storo.rep.cloaca.sample.info.tab$Color[order(Storo.rep.cloaca.sample.info.tab$Host)])) + 
  theme_bw() + theme(legend.position="top")

#Gills
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.gills.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.gills.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host != "Rhizoprionodon_terraenovae" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Sarda_sarda" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_melanopterus" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_leucas" &
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Seawater" &
                                                                                 Storo.replicate.sample.IDs.tab$Isolation_Source == "gills"]
Storo.genus.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.gills.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.gills.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.genus.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g$Color[order(Storo.genus.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Gills replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Gills replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.species.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g$Color[order(Storo.species.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Gills replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.gills.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Gills replicate samples only")

Storo.rep.gills.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.gills.replicate.sample.IDs]))
Storo.rep.gills.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.gills.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.gills.euc.dist, Storo.rep.gills.sample.info.tab$Host))
adonis2(Storo.rep.gills.euc.dist~Storo.rep.gills.sample.info.tab$Host)
Storo.rep.gills.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.gills.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.gills.sample.info.tab.phy <- sample_data(Storo.rep.gills.sample.info.tab)
Storo.rep.gills.vst.physeq <- phyloseq(Storo.rep.gills.vst.count.phy, Storo.rep.gills.sample.info.tab.phy)
Storo.rep.gills.vst.pcoa <- ordinate(Storo.rep.gills.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.gills.eigen.vals <- Storo.rep.gills.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.gills.vst.physeq, Storo.rep.gills.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.rep.gills.eigen.vals[2]/Storo.rep.gills.eigen.vals[1])) + ggtitle("PCoA - Gill replicate samples only") + 
  scale_color_manual(values=unique(Storo.rep.gills.sample.info.tab$Color[order(Storo.rep.gills.sample.info.tab$Host)])) + 
  theme_bw() + theme(legend.position="top")

#Skin
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.skin.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.skin.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host != "Rhizoprionodon_terraenovae" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Sarda_sarda" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_melanopterus" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_leucas" &
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Seawater" &
                                                                                 Storo.replicate.sample.IDs.tab$Isolation_Source == "skin"]
Storo.genus.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.skin.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.skin.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.genus.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g$Color[order(Storo.genus.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Skin replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Skin replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.species.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g$Color[order(Storo.species.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Skin replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.skin.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Skin replicate samples only")

Storo.rep.skin.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.skin.replicate.sample.IDs]))
Storo.rep.skin.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.skin.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.skin.euc.dist, Storo.rep.skin.sample.info.tab$Host))
adonis2(Storo.rep.skin.euc.dist~Storo.rep.skin.sample.info.tab$Host)
Storo.rep.skin.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.skin.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.skin.sample.info.tab.phy <- sample_data(Storo.rep.skin.sample.info.tab)
Storo.rep.skin.vst.physeq <- phyloseq(Storo.rep.skin.vst.count.phy, Storo.rep.skin.sample.info.tab.phy)
Storo.rep.skin.vst.pcoa <- ordinate(Storo.rep.skin.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.skin.eigen.vals <- Storo.rep.skin.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.skin.vst.physeq, Storo.rep.skin.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=2) + 
  coord_fixed(sqrt(Storo.rep.skin.eigen.vals[2]/Storo.rep.skin.eigen.vals[1])) + ggtitle("PCoA - Skin Samples") + 
  scale_color_manual(values=unique(Storo.rep.skin.sample.info.tab$Color[order(Storo.rep.skin.sample.info.tab$Host)])) + 
  theme_bw()

#Teeth
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.teeth.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.teeth.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host != "Rhizoprionodon_terraenovae" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Sarda_sarda" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_melanopterus" & 
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Carcharhinus_leucas" &
                                                                                 Storo.replicate.sample.IDs.tab$Host != "Seawater" &
                                                                                 Storo.replicate.sample.IDs.tab$Isolation_Source == "Teeth"]
Storo.genus.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.teeth.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.teeth.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.genus.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g$Color[order(Storo.genus.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Teeth replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Teeth replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Host)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(Storo.species.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g$Color[order(Storo.species.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g$Isolation_Source)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Teeth replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.teeth.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Teeth replicate samples only")

Storo.rep.teeth.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.teeth.replicate.sample.IDs]))
Storo.rep.teeth.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.teeth.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.teeth.euc.dist, Storo.rep.teeth.sample.info.tab$Host))
adonis2(Storo.rep.teeth.euc.dist~Storo.rep.teeth.sample.info.tab$Host)
Storo.rep.teeth.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.teeth.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.teeth.sample.info.tab.phy <- sample_data(Storo.rep.teeth.sample.info.tab)
Storo.rep.teeth.vst.physeq <- phyloseq(Storo.rep.teeth.vst.count.phy, Storo.rep.teeth.sample.info.tab.phy)
Storo.rep.teeth.vst.pcoa <- ordinate(Storo.rep.teeth.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.teeth.eigen.vals <- Storo.rep.teeth.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.teeth.vst.physeq, Storo.rep.teeth.vst.pcoa, color="Host") + 
  labs(col="Host") + geom_point(size=2) + 
  coord_fixed(sqrt(Storo.rep.teeth.eigen.vals[2]/Storo.rep.teeth.eigen.vals[1])) + ggtitle("PCoA - Teeth Samples") + 
  scale_color_manual(values=unique(Storo.rep.teeth.sample.info.tab$Color[order(Storo.rep.teeth.sample.info.tab$Host)])) + 
  theme_bw()




##HOSTS##
Storo.host.anova <- anova(betadisper(Storo.euc.dist, Storo.replicate.sample.IDs.tab$Host))

#CPe
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.CPe.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.CPe.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host == "Carcharhinus_perezii"]
Storo.genus.filt.major.taxa.proportions.CPe.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CPe.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.CPe.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CPe.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.CPe.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Carcharhinus perezi replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.CPe.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus perezi replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.CPe.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Carcharhinus perezi replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.CPe.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus perezi replicate samples only")

Storo.rep.CPe.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.CPe.replicate.sample.IDs]))
Storo.rep.CPe.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.CPe.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.CPe.euc.dist, Storo.rep.CPe.sample.info.tab$Isolation_Source))
adonis2(Storo.rep.CPe.euc.dist~Storo.rep.CPe.sample.info.tab$Isolation_Source)
Storo.rep.CPe.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.CPe.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.CPe.sample.info.tab.phy <- sample_data(Storo.rep.CPe.sample.info.tab)
Storo.rep.CPe.vst.physeq <- phyloseq(Storo.rep.CPe.vst.count.phy, Storo.rep.CPe.sample.info.tab.phy)
Storo.rep.CPe.vst.pcoa <- ordinate(Storo.rep.CPe.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.CPe.eigen.vals <- Storo.rep.CPe.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.CPe.vst.physeq, Storo.rep.CPe.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation_Source") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.rep.CPe.eigen.vals[2]/Storo.rep.CPe.eigen.vals[1])) + ggtitle("PCoA - Carcharhinus perezi replicate samples only") + 
  theme_bw() + theme(legend.position="top")

#CPl
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.CPl.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.CPl.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host == "Carcharhinus_plumbeus"]
Storo.genus.filt.major.taxa.proportions.CPl.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CPl.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.CPl.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.CPl.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.CPl.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Carcharhinus plumbeus replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.CPl.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus plumbeus replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.CPl.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Carcharhinus plumbeus replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.CPl.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Carcharhinus plumbeus replicate samples only")

Storo.rep.CPl.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.CPl.replicate.sample.IDs]))
Storo.rep.CPl.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.CPl.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.CPl.euc.dist, Storo.rep.CPl.sample.info.tab$Isolation_Source))
adonis2(Storo.rep.CPl.euc.dist~Storo.rep.CPl.sample.info.tab$Isolation_Source)
Storo.rep.CPl.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.CPl.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.CPl.sample.info.tab.phy <- sample_data(Storo.rep.CPl.sample.info.tab)
Storo.rep.CPl.vst.physeq <- phyloseq(Storo.rep.CPl.vst.count.phy, Storo.rep.CPl.sample.info.tab.phy)
Storo.rep.CPl.vst.pcoa <- ordinate(Storo.rep.CPl.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.CPl.eigen.vals <- Storo.rep.CPl.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.CPl.vst.physeq, Storo.rep.CPl.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation_Source") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.rep.CPl.eigen.vals[2]/Storo.rep.CPl.eigen.vals[1])) + ggtitle("PCoA - Carcharhinus plumbeus replicate samples only") + 
  theme_bw() + theme(legend.position="top")

#GCi
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.GCi.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.GCi.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host == "Ginglymostoma_cirratum"]
Storo.genus.filt.major.taxa.proportions.GCi.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.GCi.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.GCi.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.GCi.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.GCi.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Ginglymostoma cirratum replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.GCi.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Ginglymostoma cirratum replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.GCi.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Ginglymostoma cirratum replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.GCi.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Ginglymostoma cirratum replicate samples only")

Storo.rep.GCi.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.GCi.replicate.sample.IDs]))
Storo.rep.GCi.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.GCi.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.GCi.euc.dist, Storo.rep.GCi.sample.info.tab$Isolation_Source))
adonis2(Storo.rep.GCi.euc.dist~Storo.rep.GCi.sample.info.tab$Isolation_Source)
Storo.rep.GCi.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.GCi.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.GCi.sample.info.tab.phy <- sample_data(Storo.rep.GCi.sample.info.tab)
Storo.rep.GCi.vst.physeq <- phyloseq(Storo.rep.GCi.vst.count.phy, Storo.rep.GCi.sample.info.tab.phy)
Storo.rep.GCi.vst.pcoa <- ordinate(Storo.rep.GCi.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.GCi.eigen.vals <- Storo.rep.GCi.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.GCi.vst.physeq, Storo.rep.GCi.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.rep.GCi.eigen.vals[2]/Storo.rep.GCi.eigen.vals[1])) + ggtitle("PCoA - Skin replicate samples only") + 
  theme_bw() + theme(legend.position="top")

#GCu
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.GCu.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.GCu.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host == "Galeocerdo_cuvier"]
Storo.genus.filt.major.taxa.proportions.GCu.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.GCu.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.GCu.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.GCu.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.GCu.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Galeocerdo cuvier replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.GCu.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Galeocerdo cuvier replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.GCu.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Galeocerdo cuvier replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.GCu.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Galeocerdo cuvier replicate samples only")

Storo.rep.GCu.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.GCu.replicate.sample.IDs]))
Storo.rep.GCu.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.GCu.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.GCu.euc.dist, Storo.rep.GCu.sample.info.tab$Isolation_Source))
adonis2(Storo.rep.GCu.euc.dist~Storo.rep.GCu.sample.info.tab$Isolation_Source)
Storo.rep.GCu.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.GCu.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.GCu.sample.info.tab.phy <- sample_data(Storo.rep.GCu.sample.info.tab)
Storo.rep.GCu.vst.physeq <- phyloseq(Storo.rep.GCu.vst.count.phy, Storo.rep.GCu.sample.info.tab.phy)
Storo.rep.GCu.vst.pcoa <- ordinate(Storo.rep.GCu.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.GCu.eigen.vals <- Storo.rep.GCu.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.GCu.vst.physeq, Storo.rep.GCu.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.rep.GCu.eigen.vals[2]/Storo.rep.GCu.eigen.vals[1])) + ggtitle("PCoA - Galeocerdo cuvier replicate samples only") + 
  theme_bw() + theme(legend.position="top")

#NB
Storo.replicate.sample.IDs.tab <- data.frame(Storo.sample.info.tab)
Storo.NB.replicate.sample.IDs <- Storo.replicate.sample.IDs.tab
Storo.NB.replicate.sample.IDs <- row.names(Storo.replicate.sample.IDs.tab)[Storo.replicate.sample.IDs.tab$Host == "Negaprion_brevirostris"]
Storo.genus.filt.major.taxa.proportions.NB.replicates.only.tab.for.plot.g <- Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2[Storo.genus.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.NB.replicate.sample.IDs, ]
Storo.species.filt.major.taxa.proportions.NB.replicates.only.tab.for.plot.g <- Storo.species.filt.major.taxa.proportions.tab.for.plot.g2[Storo.species.filt.major.taxa.proportions.tab.for.plot.g2$Sample %in% Storo.NB.replicate.sample.IDs, ]
ggplot(Storo.genus.filt.major.taxa.proportions.NB.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major Genus Taxa", title="Negaprion brevirostris replicate samples only")
ggplot(Storo.genus.filt.major.taxa.proportions.NB.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Negaprion brevirostris replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.NB.replicates.only.tab.for.plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + 
  geom_jitter(aes(color=factor(Isolation_Source)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) +
  labs(x="Major ASV Taxa", title="Negaprion brevirostris replicate samples only")
ggplot(Storo.species.filt.major.taxa.proportions.NB.replicates.only.tab.for.plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", title="Negaprion brevirostris replicate samples only")

Storo.rep.NB.euc.dist <- dist(t(Storo.vst.trans.count.tab[ , colnames(Storo.vst.trans.count.tab) %in% Storo.NB.replicate.sample.IDs]))
Storo.rep.NB.sample.info.tab <- Storo.sample.info.tab[row.names(Storo.sample.info.tab) %in% Storo.NB.replicate.sample.IDs, ]
anova(betadisper(Storo.rep.NB.euc.dist, Storo.rep.NB.sample.info.tab$Isolation_Source))
adonis2(Storo.rep.NB.euc.dist~Storo.rep.NB.sample.info.tab$Isolation_Source)
Storo.rep.NB.vst.count.phy <- otu_table(Storo.vst.trans.count.tab[, colnames(Storo.vst.trans.count.tab) %in% Storo.NB.replicate.sample.IDs], taxa_are_rows=T)
Storo.rep.NB.sample.info.tab.phy <- sample_data(Storo.rep.NB.sample.info.tab)
Storo.rep.NB.vst.physeq <- phyloseq(Storo.rep.NB.vst.count.phy, Storo.rep.NB.sample.info.tab.phy)
Storo.rep.NB.vst.pcoa <- ordinate(Storo.rep.NB.vst.physeq, method="MDS", distance="euclidean")
Storo.rep.NB.eigen.vals <- Storo.rep.NB.vst.pcoa$values$Eigenvalues
plot_ordination(Storo.rep.NB.vst.physeq, Storo.rep.NB.vst.pcoa, color="Isolation_Source") + 
  labs(col="Isolation Source") + geom_point(size=1) + 
  coord_fixed(sqrt(Storo.rep.NB.eigen.vals[2]/Storo.rep.NB.eigen.vals[1])) + ggtitle("PCoA - Negaprion brevirostris replicate samples only") + 
  theme_bw() + theme(legend.position="top")
