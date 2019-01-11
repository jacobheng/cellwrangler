#'Plot mean expression heatmap
#'
#' plot_mean_exprs_heatmap() plots a heatmap showing the mean expression of a gene in a monocle
#' CellDataSet object.
#'
#' @param genes a vector of gene name(s) to plot e.g. c("Actb", "Aldoa")
#' @param cds a monocle CellDataSet object
#' @param group column name in pData(cds) to group the bar plots by on the horizontal axis.
#' @param scale_method scaling to be done for each gene; possible values are "z-score" and "rescale"; "z-score"
#' scales the values to the z-scores; "rescale" scales the values to have a minimum of 0 and maximum of 1.
#' @param cluster_groups logical; if groups should be clustered.
#' @param order_subgroups logical; if subgroups should be ordered according to clustering of first subgroup. 
#' Subgroups can be specified in a string with separation by " - ".
#' @param cluster_groups_distance distance measure used in clustering groups. Possible values are "correlation",
#' "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski" (defined in dist() function.)
#' @param cluster_groups_method clustering method used for groups. Accepts the same values as hclust: "ward.D", 
#' "ward.D2", "single", "complete", "average" (=UPGMA), "mcquitty" (=WPGMA), "median" (=WPGMC) or 
#' "centroid" (=UPGMC).
#' @param cluster_genes logical; if genes should be clustered. If TRUE, genes with zero sd will be removed.
#' @param cluster_genes_vector a vector specifying the order of genes. If NULL, genes will be clustered
#' according to distance and method specified below.
#' @param cluster_genes_distance distance measure used in clustering genes. Possible values are the
#' same as cluster_groups_distance().
#' @param cluster_genes_method clustering method used for genes. Possible values are the
#' same as cluster_groups_method().
#' @param limits limits of the color scale for the heatmap.
#' @param color_scale color scale vector for heatmap of length 3 for low, high, and na values respectively;
#' if null, limits of color scale will be determined by the squish() function
#' @param square_tiles logical; if tiles should be made square
#' 
#' @keywords plot_pattern_cor_heatmap
#' @export
#' @return a ggplot object
#' @examples
#' plot_mean_exprs_heatmap(c("Actb","Aldoa"), cds=dat)

plot_mean_exprs_heatmap <- function (genes, cds, group, scale_method = NULL, cluster_groups = F, order_subgroups = F, 
                                     cluster_groups_distance = "euclidean", cluster_groups_method = "complete", 
                                     cluster_genes = F, cluster_genes_vector = NULL, 
                                     cluster_genes_distance = "correlation", 
                                     cluster_genes_method = "complete", 
                                     color_scale = c("midnightblue", "gold", "grey50"), 
                                     limits = NULL, square_tiles = T) 
{
  #Create gene ref
  gene_ref <- as.data.frame(findGeneID(genes, cds))
  colnames(gene_ref) <- c("gene_id")
  gene_ref$gene_short_name <- findGeneName(gene_ref$gene_id, cds, unique = F)
  gene_ref$gene_short_name <- make.names(gene_ref$gene_short_name, unique = T)
  
  #subset cds
  cds_subset <- cds[findGeneID(genes, cds), ]
  #Create group_vector
  group_vector <- as.character(unique(pData(cds_subset)[, group]))
  
  if (length(genes) == 1) {
    tmp <- lapply(group_vector, function(x) {
      celltype_means <- Matrix::rowMeans(exprs(cds_subset[, 
                                                          pData(cds_subset)[, group] %in% x]))
      return(celltype_means)
    })
    tmp <- as.data.frame(matrix(unlist(tmp), nrow = 1))
    #Change rownames
    rownames(tmp) <- gene_ref$gene_short_name[match(rownames(tmp), gene_ref$gene_id)]
    #Scale data
    if(is.null(scale_method) == F) {
    if (scale_method == "z-score") {
      tmp_scale <- as.data.frame(t(scale(t(tmp))[1:length(group_vector)]))
    } else { if(scale_method == "rescale"){
      tmp_scale <- t(apply(tmp,1,FUN=rescale, to=c(0,1)))
    } else { stop("scale_method can only be z-score or rescale") }
    } } else { tmp_scale <- tmp }
    
    colnames(tmp_scale) <- group_vector
    tmp_melt <- melt(tmp_scale)
    colnames(tmp_melt) <- c("group", "Mean exprs per cell")
    tmp_melt$gene_id <- rownames(fData(cds_subset))
  } else {
    tmp <- sapply(group_vector, function(x) {
      celltype_means <- Matrix::rowMeans(exprs(cds_subset[, 
                                                          pData(cds_subset)[, group] %in% x]))
      return(celltype_means)
    })
    #Change rownames
    rownames(tmp) <- gene_ref$gene_short_name[match(rownames(tmp), gene_ref$gene_id)]
    #Scale data
    if(is.null(scale_method) == F) {
    if(scale_method == "z-score") {
      tmp_scale <- t(apply(tmp, 1, scale))
      colnames(tmp_scale) <- group_vector
    }
    else { if(scale_method == "rescale"){
      tmp_scale <- t(apply(tmp,1,FUN=rescale, to=c(0,1)))
    } else { stop("scale_method can only be z-score or rescale") }
    } } else { tmp_scale <- tmp }
    
    tmp_melt <- melt(tmp_scale)
    colnames(tmp_melt) <- c("gene", "group", "Mean exprs per cell")
  }

  tmp_melt$gene <- factor(tmp_melt$gene, levels = rev(gene_ref$gene_short_name))
  
  
  if (cluster_genes == T) {
   
    if (is.null(cluster_genes_vector) == F) {
      tmp_melt$gene <- factor(tmp_melt$gene, levels = unique(tmp_melt$gene)[cluster_genes_vector])
    }
    else {
      #Remove NAs and genes with zero sd
      tmp_scale <- tmp_scale[apply(tmp_scale, 1, sd, na.rm=TRUE) > 0,]
      #Cluster genes
      gene.order <- order.dendrogram(as.dendrogram(pheatmap:::cluster_mat((tmp_scale), 
                                                                          distance = cluster_genes_distance, method = cluster_genes_method)))
      tmp_melt$gene <- factor(tmp_melt$gene, levels = unique(tmp_melt$gene)[gene.order])
    }
  }
  else {
    tmp_melt$gene <- tmp_melt$gene
  }
  
  if (cluster_groups == T) {
  
    #Create ref_table
    ref_table <- as.data.frame(group_vector)
    colnames(ref_table) <- c("group_vector")
    #Cluster groups
    group.order <- order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(t(tmp_scale), distance = cluster_groups_distance, 
                                                                         method = cluster_groups_method)))

    ref_table$group_vector <- factor( ref_table$group_vector, 
                                      levels = ref_table$group_vector[group.order] )
    
    if (order_subgroups == T) {
      n_subgroups <- length(stringr::str_split(group_vector, 
                                               pattern = " - ")[[1]])
      for (i in 1:n_subgroups) {
        ref_table[paste("subgroup", i, sep = "")] <- as.factor(stringr::str_split_fixed(group_vector, 
                                                                                        pattern = " - ", n = n_subgroups)[, i])
      }
      for (i in (n_subgroups + 1):2) {
        ref_table <- ref_table[with(ref_table, order(ref_table[,i])), ]
      }
      group_vector_ordered <- ref_table$group_vector
      tmp_melt$group <- factor(tmp_melt$group, 
                               levels = rev(group_vector_ordered))
    }
    else {
      
      tmp_melt$group<- factor(tmp_melt$group, 
                              levels = rev(ref_table$group_vector[group.order]))
    }
  }
  else {
    tmp_melt <- tmp_melt
  }
  p <- ggplot(tmp_melt, aes(x = group, y = gene, fill = `Mean exprs per cell`))
  p <- p + geom_tile(size = 0.25) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_fill_gradient(name = "Expression", low = color_scale[1], high = color_scale[2], 
                        na.value = color_scale[3], limits = limits, oob = scales::squish)
  if (square_tiles == T) {
    p <- p + coord_equal()
  }
  else {
    p <- p
  }
  return(p)
}