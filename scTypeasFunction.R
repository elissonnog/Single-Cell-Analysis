
#working directory
setwd("/Volumes/projects_secondary/bras/elisson/projects/asap/ml_algorithms/sctype/test/")

#reading seurat object
seurat_hip <- readRDS("/Volumes/projects/data/asap_data/PD_HIP_snRNAseq_with_redos/seurat_anno.rds")
seurat_mfg <- readRDS("/Volumes/projects/data/asap_data/PD_MFG_snRNAseq_with_redos/seurat_anno.rds")
seurat_sn <- readRDS("/Volumes/projects/data/asap_data/PD_SN_snRNAseq_with_redos/seurat_anno.rds")


seurat_sn_test = runsctype(seurat_obj = seurat_sn,
                           dir = "./figures/", 
                           db = "cr", #sc(scType), kw (Kaitlyn), cr (Curated), 
                           fig.name = unique(seurat_sn@meta.data$BRAIN_REGION) #or anyname you want
)



#
###########
#Functions
#########
runsctype <- function(seurat_obj, dir, fig.name, db){
  
 
  
  #checking file 
  if (file.exists("./figures")) {
    cat("Figures Results Already Exists")
  } else {
    cat("Creating scType Figures Directory")
    dir.create("./figures")
  }
  seurat_obj$cluster_annotation = NULL
  
  # load libraries and functions
  lapply(c("dplyr","Seurat","HGNChelper","openxlsx", "patchwork","ggraph","igraph","tidyverse", "data.tree","ggplot2","hrbrthemes","viridis","tidytext"), library, character.only = T)
  
  
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  
  #DB
  if(db == "cr"){
   
    print("Using curated DB")
    db_ = ("/Volumes/projects_secondary/bras/elisson/projects/db/ScTypeDB_curated.xlsx") #curated
    
  }else if(db == "kw"){
   
    print("Using KW DB")
    db_ = ("/Volumes/projects_secondary/bras/elisson/projects/db/ScTypeDB_kwType_simple.xlsx") #kw
    #db_ = ("/Volumes/projects_secondary/bras/elisson/projects/db/ScTypeDB_kwType_full.xlsx") #kw
  }else{
  
    print("Using scType default DB")
    db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; #ORIGINAL
    
  }
  
  tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  
  # prepare gene sets
  gs_list = gene_sets_prepare(db_, tissue)
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = seurat_obj[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
  # In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
  # or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  
  
  #annotation
  seurat_obj@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat_obj@meta.data$customclassif[seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  myplot = DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') +ggtitle("scType")
  
  png(file=paste0(dir, fig.name,"_umap_annot.png"), 
      width=800, 
      height=800)
  
  print(myplot)
  dev.off()
  
  DefaultAssay(seurat_obj) <- "RNA"
  # make this the Idents()
  Idents(seurat_obj) <- seurat_obj@meta.data$customclassif
  
  
  #changing name
  seurat_obj$cluster_annotation = seurat_obj@meta.data$customclassif
  
  # prepare edges
  cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL
  
  
  # prepare nodes - if need increase the number of colours
  nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
  ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3", "#c56133", "#c3a5b4", "#37294f","#991919", "#8ad8e8", "#ffc413", "#5d4c86", "#c3a5b4", "#946aa2", "#EBA05F")
  
  
  for (i in 1:length(unique(cL_resutls$cluster))){
    dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
  }
  
  nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
  files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
  nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
  
  
  
  # Remove duplicates based on cluster column
  # test for solve duplicated data
  for( i in unique(nodes$realname)){
    a<-subset(nodes, realname == i)
    if(length(unique(a$shortName)) > 1){
      nodes[nodes["realname"] == i,"shortName"]<-min(char(unique(a$shortName)))
    }
  }
  nodes<-nodes[!duplicated(nodes),]
  
  
  mygraph <- graph_from_data_frame(edges, vertices=nodes)
  
  
  # Make the graph  
  gggr = ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
    geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
    theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
  
  
  #scater::multiplot(DimPlot(seurat_obj, reduction = "umap", label = T, repel = T, cols = ccolss), gggr, cols = 2)
  png(file=paste0(dir,fig.name,"_buble.png"), 
      width=800, 
      height=800)
  
  print(gggr)
  dev.off()
  
  #scater::multiplot(DimPlot(seurat_obj, reduction = "umap", label = T, repel = T, cols = ccolss), gggr, cols = 2)
  
  #plotting scores as bar plot
  #plot scores
  barplot = ggplot(cL_resutls %>%
           mutate(cluster = as.factor(cluster),
                  type = reorder_within(type,  n():1, cluster)),
         aes(y =type, x=scores, fill =type))+
    geom_bar(position="dodge", stat="identity")+
    #scale_x_log10() +
    scale_y_reordered() +
    ylab("Scores") + 
    xlab("Clusters")+
    theme(legend.position = "none")+
    facet_wrap(~cluster, scales = "free") #data divided and free scale by plot
  
  png(file=paste0(dir,fig.name,"_barplot.png"), 
      width=1745, 
      height=934)
  
  print(barplot)
  dev.off()
  
  return(seurat_obj)
}

