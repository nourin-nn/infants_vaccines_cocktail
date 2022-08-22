library("dplyr")
library("Seurat")
library("pheatmap")
library("ggplot2")
# ===================================
load(file.path(analysis_path, "sobj_subC_BC.RData"))
load(file.path(analysis_path, "sobj_subC_CD4TC.RData"))
load(file.path(analysis_path, "sobj_subC_CD8TC.RData"))
load(file.path(analysis_path, "sobj_subC_NKC.RData"))
load(file.path(analysis_path, "sobj_subC_MCYTE.RData"))
load(file.path(analysis_path, "sobj_subC_cDC.RData"))
load(file.path(analysis_path, "sobj_subC_pDC.RData"))
load(file.path(analysis_path, "sobj_subC_PC.RData"))
load(file.path(analysis_path, "sobj_subC_MGKC.RData"))
load(file.path(analysis_path, "sobj_subC_ERYT.RData"))
annotations <- read.csv(file = file.path(meta_data_path, "annotation.txt"))
# ===================================
sobj_ls <- list("BC" = sobj_bc,
                "CD4" = sobj_cd4tc,
                "CD8" = sobj_cd8tc,
                "NKC" = sobj_nkc,
                "Mono" = sobj_mcyte,
                "PC" = sobj_pc,
                "cDC" = sobj_cdc,
                "pDC" = sobj_pdc,
                "Mgk" = sobj_mgkc,
                "Eryt" = sobj_eryt)
subC_levels <- unlist(lapply(sobj_ls, function(x) levels(x)), use.names = FALSE)
n_subC <- length(subC_levels)
# ===================================
if (TRUE) {
  # ===================================
  genes_df_all <- data.frame()
  for (i in seq_along(sobj_ls)) {
    genes <- unique(mod$gene[grepl("IFN_", mod$ind)])

    mt_data <- FetchData(sobj_ls[[i]], c("sub_seurat_clusters", "time_point")) %>%
      tibble::rownames_to_column(var="barcode")
    
    genes_df <- FetchData(sobj_ls[[i]], genes, slot = "data") %>%
      tibble::rownames_to_column(var = "barcode") %>%
      tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
      base::merge(., mt_data, by = c("barcode"), all = TRUE) %>%
      dplyr::mutate(cpm = expm1(log1pcpm)) %>% # back from log1pCPM to CPM (cell normalized)
      dplyr::group_by(sub_seurat_clusters, time_point) %>%
      dplyr::mutate(cell_ct = length(unique(barcode))) %>% # number of cells in the cluster
      dplyr::ungroup() %>%
      dplyr::group_by(sub_seurat_clusters, time_point, gene) %>%
      dplyr::summarise(cell_ct = unique(cell_ct), 
                       cell_exp_ct = sum(cpm > 0), # number of cells with detectable (>0) expression of that gene
                       mean_cpm = mean(cpm),
                       sd_cpm = sd(cpm),
                       mean_log1pcpm = mean(log1pcpm),
                       sd_log1pcpm = sd(log1pcpm)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cell_exp_ct = ifelse(cell_exp_ct == 0, cell_exp_ct + 1, cell_exp_ct),
                    frac = cell_exp_ct / cell_ct)  
    
    genes_df_all <- genes_df_all %>%
      dplyr::bind_rows(genes_df)
  }
  table(genes_df_all$sub_seurat_clusters)
  # ===================================
  genes_df.1 <- genes_df_all %>%
    dplyr::select(sub_seurat_clusters, time_point, gene, cell_exp_ct) %>%
    tidyr::spread(key = time_point, value = cell_exp_ct, fill = 0) %>%
    dplyr::rename("count_0d" = "0d",
                  "count_7d" = "7d",
                  "count_2m" = "2m")
  
  genes_df.2 <- genes_df_all %>%
    dplyr::select(sub_seurat_clusters, time_point, gene, frac) %>%
    tidyr::spread(key = time_point, value = frac, fill = 0) %>%
    dplyr::rename("frac_0d" = "0d",
                  "frac_7d" = "7d",
                  "frac_2m" = "2m")
  
  genes_df.3 <- genes_df_all %>%
    dplyr::select(sub_seurat_clusters, time_point, gene, mean_cpm) %>%
    tidyr::spread(key = time_point, value = mean_cpm, fill = 0) %>%
    dplyr::rename("mean_cpm_0d" = "0d",
                  "mean_cpm_7d" = "7d",
                  "mean_cpm_2m" = "2m")
  
  genes_df.f <- merge(genes_df.1, genes_df.2, by = c("sub_seurat_clusters", "gene"), all = TRUE) 
  genes_df.f <- merge(genes_df.f, genes_df.3, by = c("sub_seurat_clusters", "gene"), all = TRUE) 
  genes_df.f <- genes_df.f %>%
    dplyr::mutate(score_0d = frac_0d,
                  score_7d = sqrt(frac_7d^2 + (mean_cpm_7d/mean_cpm_0d)^2), 
                  score_2m = sqrt(frac_2m^2 + (mean_cpm_2m/mean_cpm_0d)^2))

  expand_cnt <- 10
  expand_cpm <- 50
  expand_frac <- 1e-2 # filter genes with a proportion larger than 1% at B or C 
  expand_fold <- 1.5 # filter genes with a fold change larger than 2 at B or C

  genes_df.f <- genes_df.f %>%
    dplyr::mutate(expanded = (count_7d + count_2m >= expand_cnt & mean_cpm_7d + mean_cpm_2m > expand_cpm) &
                    ( (frac_7d >= expand_frac & frac_7d/frac_0d > expand_fold & mean_cpm_7d/mean_cpm_0d > expand_fold) |
                        (frac_2m >= expand_frac & frac_2m/frac_0d > expand_fold & mean_cpm_2m/mean_cpm_0d > expand_fold) ))
}
# ===================================
if (TRUE) {
  genes_df_exp <- genes_df.f %>%
    dplyr::filter(expanded == TRUE) 
  genes_per_subC <- split(genes_df_exp$gene, genes_df_exp$sub_seurat_clusters)
  genes_all <- unique(genes_df_exp$gene)
}
# ===================================
if (TRUE) {
  stackmat_df_all <- data.frame()
  labels_df_all <- data.frame()
  expanded_genes <- data.frame()
  for (subC in subC_levels) {
    genes <- genes_per_subC[[subC]]
    if (length(genes) == 0) next
    # ===================================
    db_plot <- genes_df.f %>%
      dplyr::filter(sub_seurat_clusters == subC, expanded == TRUE) %>%
      dplyr::mutate(expanded_b = (count_7d + count_2m >= expand_cnt & mean_cpm_7d + mean_cpm_2m > expand_cpm) &
                      ( (frac_7d >= expand_frac & frac_7d/frac_0d > expand_fold & mean_cpm_7d/mean_cpm_0d > expand_fold) ),
                    expanded_c = (count_7d + count_2m >= expand_cnt & mean_cpm_7d + mean_cpm_2m > expand_cpm) &
                      ( (frac_2m >= expand_frac & frac_2m/frac_0d > expand_fold & mean_cpm_2m/mean_cpm_0d > expand_fold) )) %>%
      dplyr::mutate(time_point = case_when(expanded_b == TRUE & expanded_c == FALSE ~ "7d",
                                           expanded_b == FALSE & expanded_c == TRUE ~ "2m",
                                           expanded_b == TRUE & expanded_c == TRUE ~ "72"),
                    time_point = factor(time_point, levels = c(time_points, "72"), ordered = TRUE))
    
    labels_df <- data.frame(label = c(sum(db_plot$time_point == "7d"), 
                                      sum(db_plot$time_point == "2m"),
                                      sum(db_plot$time_point == "72")),
                            time_point = c("7d", "2m", "72"),
                            x = c(21, 25, 29),
                            y = c(15.0, 15.0, 15.0),
                            fill = c(alpha("#d95f02", 0.5),
                                     alpha("#7570b3", 0.5),
                                     alpha("#ffff33", 0.5))
                            )
    
    expanded_genes <- expanded_genes %>%
      dplyr::bind_rows(db_plot %>% 
                         dplyr::mutate(time_point = as.character(time_point)))
  
    datmat_mtx <- as.matrix(db_plot[, c("frac_0d", "frac_7d", "frac_2m")])
    rownames(datmat_mtx) <- db_plot$gene
    stackmat_df <- streamplot(mtx = datmat_mtx, xvals = c(0, 7, 30), plot = F) 
    cols <- rainbow(nrow(datmat_mtx))
    stackmat_df <- stackmat_df %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(col = cols[unique(id)+1])
    # ===================================
    stackmat_df_all <- stackmat_df_all %>%
      dplyr::bind_rows(stackmat_df %>%
                         dplyr::mutate(sub_seurat_clusters = subC))
    labels_df_all <- labels_df_all %>%
      dplyr::bind_rows(labels_df %>%
                         dplyr::mutate(sub_seurat_clusters = subC))
  }
}
# ===================================
if (TRUE) {
  fct_lvl <- unique(stackmat_df_all$sub_seurat_clusters)
  fct_lvl <- fct_lvl[grepl("CD4|CD8|NKC|BC|PC|Mono|cDC|pDC|Eryt|Mgk", fct_lvl)]
  # ===================================
  db_plot <- stackmat_df_all %>%
    dplyr::group_by(sub_seurat_clusters) %>%
    dplyr::mutate(y_shift = y - 0.5*max(y)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = fct_lvl, ordered = T))
  db_labels <- labels_df_all %>%
    dplyr::mutate(y_shift = y/2) %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = fct_lvl, ordered = T))
  
  p <- ggplot(db_plot, aes(x = x, y = y_shift)) +
    getBaseTheme() +
    theme(plot.margin = unit(c(0.1, 0.1, 0.20, 0.1), "cm"),
          strip.text = element_text(size = 24, face = "bold"),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
          ) +
    coord_cartesian(ylim = c(-8.6,8.6)) +
    scale_x_continuous(breaks = c(0, 7, 30),
                       labels = time_points) +
    geom_polygon(aes(group = id, fill=col),
                 alpha=1.0, color=NA, size=0.001, show.legend=F) +
    geom_label(data = db_labels,
               aes(x = x, y = y_shift, label = label),
               color = "black", fill = db_labels$fill, fontface = "bold", size = 7.5, show.legend = F) +
    facet_wrap(. ~ sub_seurat_clusters, ncol = 10)
  # ===================================
  db_plot <- genes_df_all %>%
    dplyr::group_by(sub_seurat_clusters) %>%
    dplyr::filter(gene %in% genes_per_subC[[unique(sub_seurat_clusters)]]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = subC_levels, ordered = T))
  
  p <- ggplot(db_plot,
              aes(x=gene, y=mean_log1pcpm, color=time_point)) +
    getBaseTheme() +
    theme(plot.title = element_blank(),
          strip.text = element_text(size = 13, face = "bold"),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(size=10, face = "bold"),
          axis.text.y = element_text(size=7, face = "bold"),
          axis.ticks.y = element_blank()) +
    xlab('Gene') + ylab('Expression level') +
    geom_segment(aes(xend = gene, yend = 0), 
                 size = 0.25, color = "black") +
    geom_point(size = 2, show.legend = F) +
    coord_flip(ylim = c(0,7)) + 
    scale_color_manual(name="", values = tp_colrs) +
    facet_wrap(. ~ sub_seurat_clusters, ncol = 5, scales = "free_y")
}
# ===================================
if (TRUE) {
  genes_mtx <- genes_df_all %>%
    dplyr::filter(gene %in% genes_all) %>%
    dplyr::mutate(id = paste(sub_seurat_clusters, time_point, sep = "_")) %>%
    dplyr::select(id, gene, mean_cpm) %>%
    tidyr::spread(., key = gene, value = mean_cpm) 
  row_names <- genes_mtx$id
  genes_mtx$id <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  # ===================================
  colnames_sorted <- unique(genes_df_all$sub_seurat_clusters)
  colnames_sorted <- colnames_sorted[mixedorder(colnames_sorted)]
  colnames_sorted <- paste(rep(colnames_sorted, each = length(time_points[1:3])), time_points[1:3], sep = "_")
  genes_mtx <- genes_mtx[rownames(genes_mtx)[mixedorder(rownames(genes_mtx))], colnames_sorted]
  # ===================================
  my_col_hm <- data.frame(time_point = unlist(lapply(strsplit(colnames_sorted, split = "_"), function(x) x[3]))) 
  rownames(my_col_hm) <- colnames_sorted
  # ===================================
  cols.cor <- cor(genes_mtx, use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between samples (columns)
  rows.cor <- cor(t(genes_mtx), use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between rows (genes)
  # ===================================
  n_subC_row <- 9
  my_row_hm_1 <- hclust(as.dist(1 - rows.cor), method = "complete")
  my_row_hm_1 <- cutree(tree = my_row_hm_1, k = n_subC_row)
  my_row_hm_1 <- data.frame(Cluster = my_row_hm_1) %>%
    tibble::rownames_to_column(var = "gene")
  
  my_row_hm_2 <- sapply(rownames(genes_mtx), function(x){
    y <- mod$module[which(mod$gene == x & grepl("IFN_", mod$signiture))][1]
      stopifnot(length(y) == 1)
    return(y)
  }, USE.NAMES = FALSE)
  my_row_hm_2 <- setNames(my_row_hm_2, rownames(genes_mtx))
  my_row_hm_2 <- data.frame(module = my_row_hm_2)
  
  db <- read.csv(file = file.path(meta_data_path, "IFN_GO.txt"),
                 sep = "\t") %>%
    dplyr::rename(gene = Genes,
                  IFN_type = Group)
  
  my_row_hm_2 <- my_row_hm_2 %>%
    tibble::rownames_to_column(var="gene") %>%
    base::merge(., db, by = "gene", all.x = TRUE) %>%
    dplyr::mutate(IFN_type = case_when(is.na(IFN_type) ~ "ND",
                                       IFN_type == "Type_I" ~ "I",
                                       IFN_type == "GAMMA" ~ "II",
                                       IFN_type == "Type_I_GAMMA" ~ "I / II")) %>%
    dplyr::rename(Type = IFN_type,
                  Module = module)
  
  
  my_row_hm <- merge(my_row_hm_1, my_row_hm_2, by = c("gene"), all = TRUE) %>%
    tibble::column_to_rownames(var = "gene") %>%
    dplyr::mutate(Module = factor(Module, levels = unique(Module)[mixedorder(unique(Module))], ordered = T),
                  Type = factor(Type, levels = unique(Type)[mixedorder(unique(Type))], ordered = T))
  # ===================================
  my_colour_hm <- list(time_point = tp_colrs[1:3],
                       Module = setNames(brewer.pal(n = length(unique(my_row_hm$Module)), name = "Paired"),
                                levels(my_row_hm$Module)),
                       Type = setNames(c("#fdc086", "#386cb0", "#f0027f", "grey70"),
                                       levels(my_row_hm$Type)),
                       Cluster = setNames(brewer.pal(n = n_subC_row, name = "Set1"),
                                          1:n_subC_row))
  # ===================================
  ph <- pheatmap(genes_mtx,
                scale = "row", # This will scale the data to Z-scores (by gene) automatically
                show_colnames = T,
                show_rownames = T,
                labels_col = unlist(lapply(strsplit(colnames(genes_mtx), split = "_"), function(x) paste(x[1], x[2], sep = "_"))),
                labels_row = rownames(genes_mtx),
                color = colorpanel(1000,"steelblue","white","firebrick"),
                clustering_method = "complete",
                clustering_distance_cols = "euclidean", 
                clustering_distance_rows = as.dist(1 - rows.cor), 
                cluster_cols = T,
                cluster_rows = T,
                fontsize_row = 10,
                fontsize_col = 10,
                gaps_col = seq(3, 3*n_subC, 3),
                cutree_rows = n_subC_row,
                cutree_cols = 13,
                treeheight_col = 25,
                treeheight_row = 25,
                silent = F,
                border_color = "grey90",
                fontsize = 8,
                annotation_col = my_col_hm,
                annotation_row = my_row_hm,
                annotation_colors = my_colour_hm
  )
}
# ===================================
if (TRUE) {
  bin_mat <- matrix(0, nrow = length(genes_all), ncol = n_subC, dimnames = list(genes_all, subC_levels))
  for (x in genes_all) {
    for (y in subC_levels) {
      if (x %in% genes_per_subC[[y]]) {
        db <- genes_df.f %>%
          dplyr::filter(sub_seurat_clusters == y, gene == x, expanded == TRUE) %>%
          dplyr::mutate(expanded_b = (count_7d + count_2m >= expand_cnt & mean_cpm_7d + mean_cpm_2m > expand_cpm) &
                          ( (frac_7d >= expand_frac & frac_7d/frac_0d > expand_fold & mean_cpm_7d/mean_cpm_0d > expand_fold) ),
                        expanded_c = (count_7d + count_2m >= expand_cnt & mean_cpm_7d + mean_cpm_2m > expand_cpm) &
                          ( (frac_2m >= expand_frac & frac_2m/frac_0d > expand_fold & mean_cpm_2m/mean_cpm_0d > expand_fold) )) %>%
          dplyr::select(expanded, expanded_b, expanded_c)
        stopifnot(nrow(db) == 1)
        stopifnot(any(c(db$expanded_b, db$expanded_c)))
        
        if (db$expanded_b == TRUE){ bin_mat[x,y] <- 1}
        if (db$expanded_c == TRUE){ bin_mat[x,y] <- 2}
        if (db$expanded_b == TRUE & db$expanded_c == TRUE){ bin_mat[x,y] <- 3}
      } else {
        bin_mat[x,y] <- 0
      }
    }
  }
  x <- sort(table(genes_df_exp$gene))
  y <- sort(rowSums(bin_mat != 0))
  bin_mat <- bin_mat[, colSums(bin_mat != 0) != 0]
  bin_mat <- bin_mat[names(sort(rowSums(bin_mat != 0), decreasing = T)), names(sort(colSums(bin_mat != 0), decreasing = T))]
  p <- pheatmap(bin_mat,
                scale = "none", # This will scale the data to Z-scores (by gene) automatically
                show_colnames = T,
                show_rownames = T,
                color = c("white","#d95f02","#7570b3", "#ffff33"),
                cluster_cols = F,
                cluster_rows = F,
                fontsize_row = 10,
                fontsize_col = 10, 
                silent = F,
                legend = FALSE,
                border_color = "grey80")
}
# ===================================
