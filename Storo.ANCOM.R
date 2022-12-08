### ANCOM Analysis ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")
library(ANCOMBC)

Storo.replicate.sample.info.tab.for.ancom <- Storo.replicate.sample.IDs.tab[!Storo.replicate.sample.IDs.tab$Host == "Rhizoprionodon_terraenovae" & 
                                                                             !Storo.replicate.sample.IDs.tab$Host == "Sarda_sarda" &
                                                                             !Storo.replicate.sample.IDs.tab$Host == "Carcharhinus_melanopterus" &
                                                                             !Storo.replicate.sample.IDs.tab$Host == "Carcharhinus_leucas" &
                                                                             !Storo.replicate.sample.IDs.tab$Host == "Seawater",]
Storo.replicate.sample.info.tab.phy <- sample_data(Storo.replicate.sample.info.tab.for.ancom)
Storo.rep.count.tab <- Storo.count.tab[,-c(6,38,104,128)] #removes null samples to fix dimension is 0 error
Storo.rep.count.tab.phy <- otu_table(Storo.rep.count.tab, taxa_are_rows=T)
Storo.rep.ASV.physeq <- phyloseq(Storo.rep.count.tab.phy, Storo.tax.tab.phy, sam_data=Storo.replicate.sample.info.tab.phy)
?ancombc
Storo.asv.ancombc = ancombc(phyloseq = Storo.rep.ASV.physeq, 
                        formula = "Isolation_Source*Host",
                        p_adj_method = "holm",
                        group = "Isolation_Source",
                        lib_cut = 1000,
                        struc_zero = TRUE, 
                        neg_lb = TRUE,
                        tol = 1e-5, 
                        max_iter = 100,
                        conserve = FALSE, 
                        alpha = 0.05, 
                        global = TRUE)

Storo.ancombc.res <- Storo.asv.ancombc$res
Storo.ancombc.res.global <- Storo.asv.ancombc$res_global
Storo.ancombc.W <- Storo.ancombc.res$W
Storo.ancombc.diffabn <- Storo.ancombc.res$diff_abn
Storo.ancombc.pval <- Storo.ancombc.res$p_val
