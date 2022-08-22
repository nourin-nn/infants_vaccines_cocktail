library("dplyr")
library("Seurat")
library("pheatmap")
library("ggplot2")
library("gplots")
library("grid")
library("gridExtra")
library("cowplot")
library("ggbeeswarm")
library("ggwordcloud")
library("Nebulosa")
library("ggVennDiagram")
library("GeneOverlap")
library("slingshot")
library("ggrepel")
# ===================================
load(file.path(analysis_path, "sobj_subC_pDC.RData"))
annotations <- read.csv(file = file.path(meta_data_path, "annotation.txt"))
sobj <- sobj_pdc
sobj$time_point <- factor(sobj$time_point, levels = time_points, ordered = TRUE)
cltyp <- "pDC"
# ===================================
if (TRUE) {
  genes <- c("GZMB", "JCHAIN", "AK128525", "SERPINF1", "ITM2C", "IRF7",
             "TCF4", "BCL11A", "LILRA4", "MZB1", "DERL3", "IL3RA", "ZFAT", 
             "NRP1", "CLEC4C", "IRF8", "SLC15A4", "BLNK", "PLD4")
  
  mt_data <- FetchData(sobj, c("time_point")) %>%
    tibble::rownames_to_column(var="barcode")
  
  genes_df <- FetchData(sobj, genes, slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    dplyr::mutate(cpm = expm1(log1pcpm)) %>% # back from log1pCPM to CPM (cell normalized)
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(count = sum(log1pcpm > 0)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::mutate(gene = factor(gene, levels = unique(gene), ordered = T),
                  time_point = case_when(time_point == "0d" ~ "1",
                                         time_point == "7d" ~ "2",
                                         time_point == "2m" ~ "3"),
                  barcode_tp = paste(barcode, time_point, sep = "_")) 
  dim(genes_df)
  
  genes_mtx <- genes_df %>%
    dplyr::select(barcode_tp, gene, log1pcpm) %>%
    tidyr::spread(., key = gene, value = log1pcpm) 
  row_names <- genes_mtx$barcode_tp
  genes_mtx$barcode_tp <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  # ===================================
  colnames_sorted <- unique(genes_df$barcode_tp)
  colnames_sorted <- colnames_sorted[order(sub('.*_', '', colnames_sorted))]
  genes_mtx <- genes_mtx[, colnames_sorted]
  genes_mtx[1:10,1:10]
  # ===================================
  cols.cor <- cor(genes_mtx, use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between samples (columns)
  rows.cor <- cor(t(genes_mtx), use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between rows (genes)
  # ===================================
  my_col_hm <- FetchData(sobj, c("time_point")) %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::mutate(time_point = case_when(time_point == "0d" ~ "1",
                                         time_point == "7d" ~ "2",
                                         time_point == "2m" ~ "3"),
                  barcode = paste(barcode, time_point, sep = "_"),
                  time_point = case_when(time_point == "1" ~ "0d",
                                         time_point == "2" ~ "7d",
                                         time_point == "3" ~ "2m")) %>%
    tibble::column_to_rownames(var = "barcode")
  # ===================================
  my_colour_hm <- list(time_point = tp_colrs[1:3])
  # ===================================
  ph <- pheatmap(genes_mtx,
                 scale = "none", # This will scale the data to Z-scores (by gene) automatically
                 show_colnames = F,
                 show_rownames = T,
                 color = colorpanel(1000,"white","red"),
                 clustering_method = "complete",
                 clustering_distance_cols = as.dist(1 - cols.cor), # or use "correlation"
                 clustering_distance_rows = as.dist(1 - rows.cor), # or use "correlation"
                 cluster_cols = F,
                 cluster_rows = F,
                 fontsize_col = 10,
                 fontsize_row = 9,
                 gaps_col = cumsum(table(my_col_hm$time_point)[c("0d", "7d")]),
                 fontsize_number = 1,
                 silent = FALSE,
                 border_color = "grey60",
                 annotation_col = my_col_hm,
                 annotation_colors = my_colour_hm
  )
}
# ===================================
if (TRUE) {
  sobj_sc <- sobj
  #- IFN
  genes_ifn <- unique(mod$gene[grepl("IFN_",mod$ind)])
  length(genes_ifn)
  sobj_sc <- AddModuleScore(
    object = sobj_sc,
    assay = "RNA", 
    features = list(genes_ifn),
    name = 'interferon',
    seed = 12345,
  )
}
# ===================================
if (TRUE) {
  db_plot <- sobj_sc@meta.data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::select(barcode, time_point, interferon1) %>%
    tidyr::gather(., key = module, value = score, -barcode, -time_point) %>%
    dplyr::mutate(time_point = factor(time_point, levels = time_points, ordered = T),
                  module = case_when(module == "interferon1" ~ "Interferon"),
                  module = factor(module, levels = unique(module), ordered = T)) %>%
    dplyr::group_by(module, time_point) %>%
    dplyr::mutate(outlier = is_outlier(score)) %>%
    dplyr::ungroup()
  
  db_label <- db_plot %>%
    dplyr::filter(module == "Interferon") %>%
    dplyr::group_by(module) %>%
    dplyr::mutate(y = 0.95) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(module, time_point) %>%
    dplyr::summarise(count = n(),
                     y = unique(y)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time_point = factor(time_point, levels = time_points, ordered = T))

  error_bar <- db_plot %>%
    dplyr::filter(module == "Interferon") %>%
    dplyr::group_by(module, time_point) %>%
    dplyr::summarise(mean_sc = mean(score, na.rm = T),
                     sd_sc = sd(score, na.rm = T)) %>%
    dplyr::ungroup()
  
  p <- ggplot(data = filter(db_plot, module == "Interferon"), 
              aes(x = time_point, y = score, color = time_point, fill = time_point)) +
    getBaseTheme() +
    theme(plot.title = element_blank(),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"), #, angle = 10, hjust = 1, vjust = 1
          strip.text = element_blank(),
          legend.text=element_text(size=12)) +
    xlab("Time point") + ylab("Score") +
    coord_cartesian(ylim = c(-0.1,1)) +
    geom_errorbar(data = error_bar, 
                  aes(x = time_point, y = mean_sc,
                      ymin = mean_sc-sd_sc/2, ymax = mean_sc+sd_sc/2),
                  width=0.30, alpha=0.9, size=0.45, show.legend = F) +
    stat_summary(fun = mean, geom = "crossbar", size = 0.6, width = 0.4, show.legend = F) +
    geom_text(data = db_label,
              aes(x = time_point, y = y, color = time_point, label = count),
              size = 6.5,  fontface = "bold", show.legend = F) +
    scale_x_discrete(limits =  time_points,
                     labels = c("0d" = "A",
                                "7d" = "B",
                                "2m" = "C")) +
    scale_color_manual(name = "",
                       values = tp_colrs) +
    scale_fill_manual(name = "",
                      values = tp_colrs) +
    guides(color = guide_legend(override.aes = list(size = 0.75))) 
  plot(p)
}
# ===================================
if (TRUE) {
  mt_data <- FetchData(sobj, c("time_point")) %>%
    tibble::rownames_to_column(var="barcode")
  
  top_genes <- FetchData(sobj, genes_ifn, slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    dplyr::mutate(cpm = expm1(log1pcpm)) %>% # back from log1pCPM to CPM (cell normalized)
    dplyr::group_by(gene, time_point) %>%
    dplyr::summarise(mean_cpm = mean(cpm)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(time_point) %>%
    dplyr::top_n(n = 30, wt = mean_cpm) %>%
    dplyr::arrange(desc(mean_cpm)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(gene)

  genes_df <- FetchData(sobj, unique(top_genes), slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    dplyr::mutate(cpm = expm1(log1pcpm)) %>% # back from log1pCPM to CPM (cell normalized)
    dplyr::group_by(gene, time_point) %>%
    dplyr::summarise(mean_cpm = mean(cpm)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time_point = factor(time_point, levels = time_points, ordered = TRUE))
  
  genes_mtx <- genes_df %>%
    tidyr::spread(., key = gene, value = mean_cpm) 
  row_names <- genes_mtx$time_point
  genes_mtx$time_point <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  # ===================================
  cols.cor <- cor(genes_mtx, use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between samples (columns)
  rows.cor <- cor(t(genes_mtx), use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between rows (genes)
  # ===================================
  my_row_hm_1 <- sapply(rownames(genes_mtx), function(x){
    if (x %in% genes_ifn) return("Interferon")
    stopifnot(x %in% genes_ctx)
  }, USE.NAMES = FALSE)
  
  my_row_hm_2 <- sapply(rownames(genes_mtx), function(x){
    id <- which(grepl(pattern = x, mod$gene))
    if (x %in% genes_ifn) return(mod$module[id][as.character(mod$signiture[id]) == "IFN_Response"][1])
  }, USE.NAMES = FALSE)
  
  my_row_hm <- data.frame(signiture = my_row_hm_1,
                          Module = my_row_hm_2)
  rownames(my_row_hm) <- rownames(genes_mtx)
  
  db <- read.csv(file = file.path(meta_data_path, "IFN_GO.txt"),
                 sep = "\t") %>%
    dplyr::rename(gene = Genes,
                  IFN_type = Group)
  
  my_row_hm <- my_row_hm %>%
    tibble::rownames_to_column(var = "gene") %>%
    base::merge(., db, by = "gene", all.x = TRUE) %>%
    dplyr::mutate(IFN_type = case_when(is.na(IFN_type) ~ "ND",
                                       IFN_type == "Type_I" ~ "I",
                                       IFN_type == "GAMMA" ~ "II",
                                       IFN_type == "Type_I_GAMMA" ~ "I / II")) %>%
    dplyr::rename(Type = IFN_type) %>%
    tibble::column_to_rownames(var = "gene") %>%
    dplyr::mutate(signiture = factor(signiture, levels = unique(signiture)[mixedorder(unique(signiture))], ordered = T),
                  Module = factor(Module, levels = unique(Module)[mixedorder(unique(Module))], ordered = T),
                  Type = factor(Type, levels = unique(Type)[mixedorder(unique(Type))], ordered = T))
  
  my_row_hm$signiture <- NULL
  # ===================================
  my_colour_hm <- list(timePoint = tp_colrs[1:3],
                       Module = setNames(brewer.pal(n = length(unique(my_row_hm$Module)), name = "Paired"),
                                         levels(my_row_hm$Module)),
                       Type = setNames(c("#fdc086", "#386cb0", "#f0027f", "grey70"),
                                       levels(my_row_hm$Type))
  )
  # ===================================
  n_clust <- 3
  # ===================================
  colnames(genes_mtx) <- c("A", "B", "C")
  p <- pheatmap(genes_mtx,
                scale = "row", # This will scale the data to Z-scores (by gene) automatically
                show_colnames = T,
                show_rownames = T,
                labels_row = rownames(genes_mtx),
                color = colorpanel(1000,"blue","white","red"),
                clustering_method = "complete",
                clustering_distance_cols = as.dist(1 - cols.cor), # or use "correlation"
                clustering_distance_rows = as.dist(1 - rows.cor), # or use "correlation"
                cluster_cols = F,
                cluster_rows = TRUE,
                fontsize_col = 12,
                fontsize_row = 14,
                cellwidth = 15, 
                cellheight = 15,
                cutree_rows = n_clust,
                cutree_cols = n_clust,
                treeheight_row = 0,
                treeheight_col = 20.0,
                fontsize_number = 1,
                silent = FALSE,
                border_color = NA,
                annotation_row = my_row_hm,
                annotation_colors = my_colour_hm,
                legend = TRUE,
                annotation_legend = TRUE
  )
  # ===================================
}

  