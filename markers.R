brain_gene <- list(
  "Oligodendrocytes" = c("Cldn11", "Cnp39","PDGFRA", "CSPG4"),
  "Oligodendrocyte progenitor cells" = c("Opc", "Pdgfr", "Olig2"),
  
  "Microlia" = c("P2ry12", "C1qa", "Cx3cr1", "A2m"),
  "Astrocytes" = c("Scl6a11", "Ntsr2"),
  
  "Neurons" = c("Gad1","Gad2","Gja1","Slc17a7","Satb2","Lmo7","Pcp4","Pcdh8"),
  "neuron" = str_to_title(str_to_lower    ( c("SLC17A7","SATB2","LMO7","PCP4","PCDH8"))),
  "Excitatory Neurons" = c("Syt1", "Snap25"),
  "Neuronal precursors" = c("Sox11", "Stmn2"),
  "Endothelial cells" = c("Itm2a49", "Flt1"),
  "GABAergic neuron" = c("Npy51,Nr2f2"),
  "Excitatory neurons" = c("Gad1", "Gad2"),
  "inhibitory neurons" = c("Gja1"),
  
  "VSMCs" = c("Myh11","Acta2","Tagln"),
  "Fibroblasts" = c("Pclaf","Col1a1","Dcn"),
  "Vascular cells" = c("Myl9", "Mgp46"),
  "Endothelial vascular cells" = c("Cldn5","Pecam1","Cdh5","Mbp","Mobp","Mog"),
  "Choroid cells" = c("Lcn2")  # 脉络膜
)

immune_gene <- list(
  "T cells" = c("Cd3d", "Cd3e", "Cd3g"),
  
  
  "B cells" = c("Cd79a","Ms4a1","Cd19"),
  
  "DCs" = c("Cd209a", "Flt3", "Klrd1"),
  
  "Macrophages" = c("Lyz2", "Cd68", "Itgam"),
  
  "Neutrophils" = c("Cxcr4", "S100a9", "S100a8")
  
)

micro_polar <- c("Nos2","Il1b","Tnf","Ccl3","Ccl2","Arg1","Ch    il3","Il10")
micro_meta <- c("Hk2","Pkm","Pfkm","Ogdh","Atp5b")
marker_genes <- list(
  "Microglia" = c("A2m","P2ry12","Itgam",
                  "Cx3cr1","Tmem119","Hexb",
                  "Sparc",
                  "Ccl3","Ccl4","Fosb","Atf3" ,"Trem2","Klf2"),
  "DC"=c("H2-Eb1","H2-Ab1","H2-Aa","Cd74","Vim"),
  "Monocytes"=c("Vim","Ly6c2","Plac8","Ifitm3","S100a4"),
  "T cell"=c("Cd3d","Trbc2","Cd3e","Trac","Ms4a1","Ms4a4b","C    cl5"),
  "B cell"=c("Ms4a1","Igkc","Cd79a","Cd79b"),
  "Neutrophil"=c("S100a8","Retnlg","S100a9","Mmp8","Mmp9"),
  "Granulo"=c("Lyz2","Fn1","Hp","S100a8"),
  "Macrophage"=c("Ifit3","Isg15","Cd74","H2-Ab1"),
  "NK"=c("Nkg7","AW112010","Klrk1","Prf1","Gzma"),
  "Neuron" = c("Gad1","Gad2","Gja1",str_to_title(str_to_l( c("SLC17A7","SATB2","LMO7","PCP4","PCDH8")))))