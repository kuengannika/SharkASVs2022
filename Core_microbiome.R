install.packages("remotes")
remotes::install_github("Russel88/MicEco", force=TRUE)

library("MicEco")
library("VennDiagram")
library(gplots)
sample_data(Storo.ASV.physeq)


skin_rel = microbiome::transform(Storo.skin.vst.physeq, "compositional")
teeth_rel = microbiome::transform(Storo.teeth.vst.physeq, "compositional")

skin_uni <- unique(as.character(meta(skin_rel)$Host))
teeth_uni <-unique(as.character(meta(teeth_rel)$Host))
skin_uni_top = skin_uni[c(1, 2, 3, 4, 7)]
teeth_uni_top = teeth_uni[c(1, 2, 4, 5, 7)]

skin_list_core <- c() 
teeth_list_core <- c()

for (n in skin_uni_top){ 
  
  ps.sub <- subset_samples(skin_rel, Host == n) 
  
  core_m <- core_members(ps.sub,  
                         detection = 0.001,  
                         prevalence = 0.75) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) 
  skin_list_core[[n]] <- core_m 
}

for (n in teeth_uni_top){ 
  
  ps.sub <- subset_samples(teeth_rel, Host == n) 
  
  core_m <- core_members(ps.sub,
                         detection = 0.001, 
                         prevalence = 0.75) 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) 
  teeth_list_core[[n]] <- core_m 
}

plot(venn(skin_list_core))
plot(venn(teeth_list_core))

library(VennDiagram)
venn.diagram(skin_list_core, filename = "skin_list_core.png")
venn.diagram(teeth_list_core)

ggVennDiagram(skin_list_core, label_alpha = 0, label = "count", set_size = 2)
ggVennDiagram(teeth_list_core, label_alpha = 0, label = "count", set_size = 3)


skin_list_core
teeth_list_core