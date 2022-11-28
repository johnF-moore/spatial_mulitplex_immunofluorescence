library(dplyr)
library(ggplot2)
library(here)
library(Biobase)
library(flowCore)
library(flowViz)
library(Seurat)
library(SeuratObject)
library(readr)
library(DT)

##--------
# Functions
##--------
write_df_to_fcs= function(df, output_path){
  column_names = dimnames(df)[[2]]
  meta <- data.frame(name= column_names, desc= column_names)
  meta$range    <- apply(apply(df[,column_names],2,range),2,diff)
  meta$minRange <- apply(df[,column_names],2,min)
  meta$maxRange <- apply(df[,column_names],2,max)
  head(meta)
  # all these are required for the following steps to work
  # a flowFrame is the internal representation of a FCS file
  matrix_exprs_data <- data.matrix(df[,column_names])
  ff <- new("flowFrame", exprs=matrix_exprs_data,
            parameters=AnnotatedDataFrame(meta))
  ## generates an object of class flowFrame w/ the given data
  # now you can save it back to the filesystem
  
  write.FCS(ff, output_path)  
}


fcs_to_df= function(fcs_path){
  df <- as.data.frame(exprs(read.FCS(fcs_path)))
  return(df)
}


cluster_info= function(df){
  props= df %>% 
    count(cluster_id) %>% 
    mutate(prop= round(n/sum(n), digits= 3),
           n_clusters= length(unique(cluster_id))) %>% 
    arrange(desc(n))
  props$cumprop= cumsum(props$prop)
  return(props)
}


fcs.df_to_seurat= function(df, col_index, scale_factor, assay= "CODEX"){
  ## Takes data.frame output from fcs_to_df() and makes a Seurat object 
  rownames(df)= paste0("bar.", rownames(df))
  df_md= df[,(col_index+1):ncol(df), drop= FALSE]
  seurat_data <- as.matrix(df[,1:col_index])
  
  ## Normalization and scaling
  seurat_data <- scale(seurat_data, 
                       center= FALSE,
                       scale = apply(seurat_data, 2, sd, na.rm = TRUE))
  seurat_data <- asinh(seurat_data/scale_factor)


  
  obj= CreateSeuratObject(counts    = seurat_data,
                          assay     = assay,
                          meta.data = df_md)
  obj$cluster_id <- as.factor(df$cluster_id)
  obj$barcode    <- rownames(df)
  
  obj <- RunUMAP(obj, 
                 features= rownames(obj@assays[["CODEX"]]@counts),
                 dims= NULL,
                 verbose= TRUE)
  return(obj)
}


seurat_figs_wrapper= function(obj, nClusters){
  ## Wrapper function to make the preliminary figures that I use to assess the success of the spatial multiplex clustering 
  fig.list= vector(mode= "list", length= 8)
  names(fig.list)= c("obj", "UMAP","facet_UMAP", "feat_UMAPs",
                     "dotplot","cluster_data","spatial_clusters","facet_spatial")
  all_feats= rownames(obj@assays[["CODEX"]]@counts)
  cluster_data= cluster_info(obj@meta.data)
  if(is.numeric(nClusters)){
    kept_clusters= cluster_data$cluster_id[1:nClusters]
  } else if(is.character(nClusters) || is.factor(nClusters)){
    kept_clusters= cluster_data$cluster_id[cluster_data$cluster_id %in% nClusters]
  }
  else{
    stop("nClusters must be of class numeric, character, or factor.")
  }
  sub_obj= subset(obj, cluster_id %in% kept_clusters)
  my_colors= colorRampPalette(brewer.pal(n= 8, "Dark2"))(length(unique(kept_clusters)))
  names(my_colors)= kept_clusters
  my_colors["167"] <- "purple"
  my_colors["86"] <- "cyan"
  my_colors["110"] <- "white"
    ## Adding specific colors for some cell types 
    ## Neither of these colors are in the palette already
    ##  any(col2hcl(col= c("cyan","purple")) %in% my_colors)
  
  fig.list[[1]] <- obj
  fig.list[[2]] <- DimPlot(sub_obj, group.by= "cluster_id", cols = my_colors)

  fig.list[[3]] <- DimPlot(sub_obj, split.by= "cluster_id", ncol= 3, cols= my_colors)
  fig.list[[4]] <- FeaturePlot(sub_obj, features= all_feats, ncol= 3) 
  Idents(sub_obj) <- "cluster_id"
  fig.list[[5]] <- DotPlot(sub_obj, features = all_feats, cols = "PiYG") +
    RotatedAxis()
  fig.list[[6]] <- cluster_info(sub_obj@meta.data)
  ## I can't get heatmap working. 
  # DoHeatmap(sub_obj, features= all_feats, group.by= "cluster_id", slot= "counts")
  
  spatial_plot <- ggplot(sub_obj@meta.data) +
    aes(x= X_centroid, y= Y_centroid, color= cluster_id) + 
    geom_point(size= 0.1) + 
    theme_minimal() + 
    xlab("") + ylab("") +
    theme(axis.text = element_blank(),
          legend.position= "top") + 
    scale_color_manual(values= my_colors) +
    coord_equal()
  fig.list[[7]] <- spatial_plot
  fig.list[[8]] <- spatial_plot + 
    facet_wrap(~cluster_id) + 
    geom_point(size= 0.5)
  
  return(fig.list)
}

fcs_seurat_pipeline= function(path=NULL, df= NULL, col_index, nClusters, downsample,
                              scale_factor= 5){
  ## This function goes from data.frame to fcs to seurat to preliminary figures 
  if(is.null(df) && is.null(path)){
    stop("df and path cannot both be NULL. One must be set")
  } else if(!is.null(df) && !is.null(path)){
    stop("df and path cannot both be set. One must be NULL.")
  } else if(is.null(df) && is.character(path)){
    df= fcs_to_df(fcs_path)  
  } else if(is.null(path) && !is.data.frame(df)){
    stop("path is NULL, but df is not a data.frame.")
  }
  df <- as.data.frame(df)
  if(downsample){ df= df %>% slice_sample(prop= 0.1) }
  
  obj = fcs.df_to_seurat(df, col_index= col_index, scale_factor= scale_factor)
  obj_figs = seurat_figs_wrapper(obj, nClusters= nClusters)
  
  return(obj_figs)
}
