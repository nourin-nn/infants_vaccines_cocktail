library("dplyr")
library("Seurat")
library("pheatmap")
library("ggplot2")
library("gplots")
library("grid")
library("gridExtra")
library("cowplot")
# ===================================
sobj_res <- readRDS(file.path(analysis_path, "final_object.rds"))
sobj <- sobj_res
# ===================================
Idents(sobj) <- "cell_type"
# ===================================
ab_response_ls <- list("Hepatitis-B" = setNames(c("FR", "FR", "NI", "FR", "FR", "FR"),
                                                c("I1", "I2", "I3", "I4", "I5", "I6")),
                       "Diptheria" = setNames(c("NR", "NR", "NI", "NR", "NR", "FR"),
                                              c("I1", "I2", "I3", "I4", "I5", "I6")),
                       "Tetanus" = setNames(c("NR", "NR", "NI", "NR", "NR", "NR"),
                                            c("I1", "I2", "I3", "I4", "I5", "I6")),
                       "Pneumococcus" = setNames(c("PR", "PR", "NI", "NR", "PR", "FR"),
                                                 c("I1", "I2", "I3", "I4", "I5", "I6")),
                       "Rotavirus" = setNames(c("FR", "FR", "NI", "FR", "NR", "FR"),
                                              c("I1", "I2", "I3", "I4", "I5", "I6")))
ab_response_df <- as.data.frame(ab_response_ls)
# ===================================
if (TRUE) {
  genes_mtx <- sobj@meta.data %>%
    dplyr::group_by(sub_seurat_clusters, time_point, subject) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = paste(subject, time_point, sep = "_")) %>%
    dplyr::select(id, sub_seurat_clusters, count) %>%
    tidyr::spread(., key = sub_seurat_clusters, value = count, fill = 0)
  row_names <- genes_mtx$id
  genes_mtx$id <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  genes_mtx <- t(t(genes_mtx)/colSums(genes_mtx)) # normalize to one fro given infants and time points
  # ===================================
  colnames_sorted <- unique(sobj@meta.data$subject)
  colnames_sorted <- colnames_sorted[mixedorder(colnames_sorted)]
  colnames_sorted <- paste(rep(colnames_sorted, each = length(time_points[1:3])), time_points[1:3], sep = "_")
  genes_mtx <- genes_mtx[rownames(genes_mtx)[mixedorder(rownames(genes_mtx))], colnames_sorted]
  row_sorted <- c(rownames(genes_mtx)[grepl("CD4", rownames(genes_mtx))], 
                  rownames(genes_mtx)[grepl("CD8", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("BC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("Mono", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("NKC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("Eryt", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("pDC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("cDC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("PC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("Mgk", rownames(genes_mtx))])
  genes_mtx <- genes_mtx[row_sorted, ]
  # ===================================
  my_col_hm <- data.frame(infant = unlist(lapply(strsplit(colnames_sorted, split = "_"), function(x) x[1])), 
                          time_point = unlist(lapply(strsplit(colnames_sorted, split = "_"), function(x) x[2])))%>%
    base::merge(., ab_response_df %>%
                  tibble::rownames_to_column(var = "infant"),
                by = "infant", all = TRUE)
  my_col_hm$infant <- NULL
  rownames(my_col_hm) <- colnames_sorted
  # ===================================
  cols.cor <- cor(genes_mtx, use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between samples (columns)
  rows.cor <- cor(t(genes_mtx), use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between rows (genes)
  # ===================================
  my_row_hm <- data.frame(Culster = rownames(genes_mtx)) 
  my_row_hm$Group <- sapply(my_row_hm$Culster, function(x) { 
    if (grepl("BC|PC|CD4|CD8|NKC", x))  return("Adaptive immune")
    if (grepl("Mono|DC", x))  return("Innate immune")
    return("Non immune")
  })
  rownames(my_row_hm) <- my_row_hm$Culster
  my_row_hm$Culster <- NULL
  # ===================================
  my_colour_hm <- list(time_point = tp_colrs[1:3],
                       Hepatitis.B = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                              c("FR", "PR", "NR", "NI")),
                       Diptheria = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                            c("FR", "PR", "NR", "NI")),
                       Tetanus = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                          c("FR", "PR", "NR", "NI")),
                       Pneumococcus = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                               c("FR", "PR", "NR", "NI")),
                       Rotavirus = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                          c("FR", "PR", "NR", "NI")),
                       Group = setNames(c("#6a51a3", "#bcbddc", "#efedf5"),
                                        c("Innate immune", "Adaptive immune", "Non immune")))
  # ===================================
  ph <- pheatmap(genes_mtx,
                 scale = "row", # This will scale the data to Z-scores (by gene) automatically
                 show_colnames = T,
                 show_rownames = T,
                 labels_col = unlist(lapply(strsplit(colnames(genes_mtx), split = "_"), function(x) x[1])),
                 labels_row = rownames(genes_mtx),
                 color = colorpanel(1000,"steelblue","white","firebrick"),
                 clustering_method = "complete",
                 clustering_distance_cols = "euclidean", 
                 clustering_distance_rows = "euclidean", 
                 cluster_cols = F,
                 cluster_rows = F,
                 fontsize_row = 10,
                 fontsize_col = 10, 
                 cellwidth = 15,
                 cellheight = 15,
                 gaps_col = seq(3, 3*5, 3),
                 gaps_row = c(5,8,13,15,17,18,19,20,21),
                 cutree_rows = n_subC_row,
                 cutree_cols = 13,
                 treeheight_col = 25,
                 treeheight_row = 25,
                 silent = F,
                 border_color = NA, 
                 fontsize = 8,
                 annotation_col = my_col_hm,
                 annotation_row = my_row_hm,
                 annotation_colors = my_colour_hm
  )
}
# ===================================
#- IFN
genes_ifn <- unique(mod$gene[grepl("IFN_",mod$ind)])
length(genes_ifn)
sobj <- AddModuleScore(
  object = sobj,
  assay = "RNA", 
  features = list(genes_ifn),
  name = 'interferon'
)
# ===================================
if (TRUE) {
  genes_mtx <- sobj@meta.data %>%
    dplyr::group_by(sub_seurat_clusters, time_point, subject) %>%
    dplyr::summarise(mean_ifn = mean(interferon1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = paste(subject, time_point, sep = "_")) %>%
    dplyr::select(id, sub_seurat_clusters, mean_ifn) %>%
    tidyr::spread(., key = sub_seurat_clusters, value = mean_ifn, fill = 0)
  row_names <- genes_mtx$id
  genes_mtx$id <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  # ===================================
  colnames_sorted <- unique(sobj@meta.data$subject)
  colnames_sorted <- colnames_sorted[mixedorder(colnames_sorted)]
  colnames_sorted <- paste(rep(colnames_sorted, each = length(time_points[1:3])), time_points[1:3], sep = "_")
  genes_mtx <- genes_mtx[rownames(genes_mtx)[mixedorder(rownames(genes_mtx))], colnames_sorted]
  row_sorted <- c(rownames(genes_mtx)[grepl("CD4", rownames(genes_mtx))], 
                  rownames(genes_mtx)[grepl("CD8", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("BC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("Mono", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("NKC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("Eryt", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("pDC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("cDC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("PC", rownames(genes_mtx))],
                  rownames(genes_mtx)[grepl("Mgk", rownames(genes_mtx))])
  genes_mtx <- genes_mtx[row_sorted, ]
  # ===================================
  my_col_hm <- data.frame(infant = unlist(lapply(strsplit(colnames_sorted, split = "_"), function(x) x[1])), 
                          time_point = unlist(lapply(strsplit(colnames_sorted, split = "_"), function(x) x[2])))%>%
    base::merge(., ab_response_df %>%
                  tibble::rownames_to_column(var = "infant"),
                by = "infant", all = TRUE)
  my_col_hm$infant <- NULL
  rownames(my_col_hm) <- colnames_sorted
  # ===================================
  cols.cor <- cor(genes_mtx, use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between samples (columns)
  rows.cor <- cor(t(genes_mtx), use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between rows (genes)
  # ===================================
  my_row_hm <- data.frame(Culster = rownames(genes_mtx)) 
  my_row_hm$Group <- sapply(my_row_hm$Culster, function(x) { 
    if (grepl("BC|PC|CD4|CD8|NKC", x))  return("Adaptive immune")
    if (grepl("Mono|DC", x))  return("Innate immune")
    return("Non immune")
  })
  rownames(my_row_hm) <- my_row_hm$Culster
  my_row_hm$Culster <- NULL
  # ===================================
  my_colour_hm <- list(time_point = tp_colrs[1:3],
                       Hepatitis.B = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                              c("FR", "PR", "NR", "NI")),
                       Diptheria = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                            c("FR", "PR", "NR", "NI")),
                       Tetanus = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                          c("FR", "PR", "NR", "NI")),
                       Pneumococcus = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                               c("FR", "PR", "NR", "NI")),
                       Rotavirus = setNames(c("#006d2c", "#66c2a4", "#edf8fb", "white"),
                                          c("FR", "PR", "NR", "NI")),
                       Group = setNames(c("#6a51a3", "#bcbddc", "#efedf5"),
                                        c("Innate immune", "Adaptive immune", "Non immune")))
  # ===================================
  ph <- pheatmap(genes_mtx,
                 scale = "row", # This will scale the data to Z-scores (by gene) automatically
                 show_colnames = T,
                 show_rownames = T,
                 labels_col = unlist(lapply(strsplit(colnames(genes_mtx), split = "_"), function(x) x[1])),
                 labels_row = rownames(genes_mtx),
                 color = colorpanel(1000,"steelblue","white","firebrick"),
                 clustering_method = "complete",
                 clustering_distance_cols = "euclidean", 
                 clustering_distance_rows = "euclidean", 
                 cluster_cols = F,
                 cluster_rows = F,
                 fontsize_row = 10,
                 fontsize_col = 10, 
                 cellwidth = 15,
                 cellheight = 15,
                 gaps_col = seq(3, 3*5, 3),
                 gaps_row = c(5,8,13,15,17,18,19,20,21),
                 cutree_rows = n_subC_row,
                 cutree_cols = 13,
                 treeheight_col = 25,
                 treeheight_row = 25,
                 silent = F,
                 border_color = NA,
                 fontsize = 8,
                 annotation_col = my_col_hm,
                 annotation_row = my_row_hm,
                 annotation_colors = my_colour_hm
  )
}
# ===================================