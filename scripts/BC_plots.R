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
library("ggrepel")
# ===================================
load(file.path(analysis_path, "sobj_subC_BC.RData"))
annotations <- read.csv(file = file.path(meta_data_path, "annotation.txt"))
sobj <- sobj_bc
sds <- sds_bc
cltyp <- "BC"
subC_colrs <- brewer.pal(n = length(levels(sobj)), name = "Set1")
subC_levels <- levels(sobj)
n_subC <- length(subC_levels)
# ===================================
if (TRUE) {
  tbl <- table(sobj@meta.data$sub_seurat_clusters, sobj@meta.data$time_point)
  db_plot <- data.frame(sub_seurat_clusters = rownames(tbl)[row(tbl)], 
                        time_point = colnames(tbl)[col(tbl)],
                        count = c(tbl)) %>%
    dplyr::group_by(sub_seurat_clusters) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(data.frame(sub_seurat_clusters = "All",
                                time_point = colnames(tbl),
                                count = colSums(tbl),
                                frac = colSums(tbl)/sum(tbl))) %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = c(subC_levels, "All"), ordered = T),
                  time_point = factor(time_point, levels = time_points, ordered = T))
  
  db_label_1 <- data.frame(sub_seurat_clusters = rownames(tbl), 
                           y = 1.0475, 
                           count = rowSums(tbl)) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::bind_rows(data.frame(sub_seurat_clusters = "All",
                                y = 1.0475,
                                count = sum(tbl),
                                frac = 1)) %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = c(subC_levels, "All"), ordered = T))
  # ===================================
  p <- ggplot() +
    getBaseTheme() +
    theme(axis.text = element_text(size = 22, face = "bold"),
          axis.title = element_text(size = 22, face = "bold"),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 20),
          legend.key.width=grid::unit(30, "points"),
          legend.spacing = grid::unit(30, "points")) +
    xlab("") + ylab("Proportion of cells") + 
    scale_y_continuous(expand = expansion(mult = c(0.01, .1))) +
    geom_bar(data = db_plot,
             aes(x = sub_seurat_clusters, y = frac, fill = time_point),
             width = 0.95, stat = "identity", position = position_stack(reverse = TRUE), show.legend = F) +
    geom_text(data = db_plot,
              aes(x = sub_seurat_clusters, y = frac, fill = time_point, label = paste0(round(frac, 3)*100, "%")),
              size = 7, fontface = "bold", position = position_stack(vjust = 0.5, reverse = TRUE)) +
    geom_text(data = db_label_1, 
              aes(x = sub_seurat_clusters, y = y, label = paste0(scales::comma(count, accuracy = 1), 
                                                                 "\n(", round(frac, 2)*100, "%)")),
              size = 5.25, hjust = 0.5, fontface = "bold", color = "firebrick", angle = 0) +
    scale_fill_manual(name = "",
                      values = tp_colrs) +
    coord_flip(ylim = c(0,1)) +
    guides(fill = guide_legend(override.aes = list(size = 8)))
  plot(p)
  # ===================================
  tbl <- table(sobj@meta.data$sub_seurat_clusters, sobj@meta.data$time_point)
  db_plot <- data.frame(sub_seurat_clusters = rownames(tbl)[row(tbl)], 
                        time_point = colnames(tbl)[col(tbl)],
                        count = c(tbl)) %>%
    dplyr::group_by(time_point) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(data.frame(sub_seurat_clusters = rownames(tbl),
                                time_point = "All",
                                count = rowSums(tbl),
                                frac = rowSums(tbl)/sum(tbl))) %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = c(subC_levels, "All"), ordered = T),
                  time_point = factor(time_point, levels = c(time_points, "All"), ordered = T))
  
  db_label_1 <- data.frame(time_point = colnames(tbl), 
                           y = 1.0475, 
                           count = colSums(tbl)) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::bind_rows(data.frame(time_point = "All",
                                y = 1.0475,
                                count = sum(tbl),
                                frac = 1)) %>%
    dplyr::mutate(time_point = factor(time_point, levels = c(time_points, "All"), ordered = T))
  # ===================================
  p <- ggplot() +
    getBaseTheme() +
    theme(axis.text = element_text(size = 30, face = "bold"),
          axis.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.key.width=grid::unit(30, "points"),
          legend.spacing = grid::unit(30, "points")) +
    xlab("") + ylab("Proportion of cells") + 
    scale_y_continuous(expand = expansion(mult = c(0.01, .1))) +
    scale_x_discrete(limits =  c(time_points, "All"),
                     labels = c("0d" = "A",
                                "7d" = "B",
                                "2m" = "C",
                                "All" = "All")) +
    geom_bar(data = db_plot,
             aes(x = time_point, y = frac, fill = sub_seurat_clusters),
             width = 0.95, stat = "identity", position = position_stack(reverse = TRUE), show.legend = F) +
    geom_text(data = filter(db_plot, frac > 0.05),
              aes(x = time_point, y = frac, fill = sub_seurat_clusters, label = paste0(round(frac, 3)*100, "%")),
              size = 7, fontface = "bold", position = position_stack(vjust = 0.5, reverse = TRUE)) +
    scale_size(range = c(3,6)) +
    guides(size = FALSE) +
    geom_text(data = db_label_1, 
              aes(x = time_point, y = y, label = paste0(scales::comma(count, accuracy = 1), 
                                                        "\n(", round(frac, 2)*100, "%)")),
              size = 5.25, hjust = 0.5, fontface = "bold", color = "firebrick", angle = 0) +
    scale_fill_manual(name = "",
                      values = setNames(subC_colrs, levels(sobj))) +
    coord_flip(ylim = c(0,1)) +
    guides(fill = guide_legend(override.aes = list(size = 6), nrow = 1))
  plot(p)
}
# ===================================
if (TRUE) {
  db_plot <- sobj@meta.data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::group_by(sample, sub_seurat_clusters) %>%
    dplyr::summarise(cell_cnt = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(., sample = unique(sample), 
                    sub_seurat_clusters = unique(sub_seurat_clusters),
                    fill = list(cell_cnt = 0)) %>%
    base::merge(., demog_tbl, by = c("sample"), all.x = TRUE) 
  
  db <- db_plot %>%
    dplyr::select(subject, time_point, sub_seurat_clusters, cell_cnt) %>%
    dplyr::mutate(subject = as.character(subject)) %>%
    tidyr::spread(., key = time_point, value = cell_cnt) %>%
    dplyr::rename(A = `0d`,
                  B = `7d`,
                  C = `2m`) %>%
    dplyr::mutate(A = ifelse(is.na(A), 0, A), 
                  B = ifelse(is.na(B), 0, B), 
                  C = ifelse(is.na(C), 0, C), 
                  # =================================
                  sudo_A = A + 1,
                  sudo_B = B + 1,
                  sudo_C = C + 1) %>%
    dplyr::group_by(subject) %>%
    dplyr::mutate(total_A = sum(A),
                  total_B = sum(B),
                  total_C = sum(C),
                  total_sudo_A = sum(sudo_A),
                  total_sudo_B = sum(sudo_B),
                  total_sudo_C = sum(sudo_C)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(frac_A = A / total_A,
                  frac_B = B / total_B,
                  frac_C = C / total_C,
                  # =================================
                  frac_sudo_A = sudo_A / total_sudo_A,
                  frac_sudo_B = sudo_B / total_sudo_B,
                  frac_sudo_C = sudo_C / total_sudo_C) %>%
    dplyr::mutate(rel_diff_B = case_when(A == 0 & B == 0 ~ as.numeric(NA),
                                         A != 0 & B == 0 ~ (frac_B - frac_A)/frac_A,
                                         A == 0 & B != 0 ~ (frac_sudo_B - frac_sudo_A)/frac_sudo_A,
                                         A != 0 & B != 0 ~ (frac_B - frac_A)/frac_A),
                  # =================================
                  rel_diff_C = case_when(A == 0 & C == 0 ~ as.numeric(NA),
                                         A != 0 & C == 0 ~ (frac_C - frac_A)/frac_A,
                                         A == 0 & C != 0 ~ (frac_sudo_C - frac_sudo_A)/frac_sudo_A,
                                         A != 0 & C != 0 ~ (frac_C - frac_A)/frac_A))
  # find significants
  sig.df <- sapply(unique(db$sub_seurat_clusters), function(x) {
    a <- db$frac_sudo_A[db$sub_seurat_clusters == x]
    b <- db$frac_sudo_B[db$sub_seurat_clusters == x]
    c <- db$frac_sudo_C[db$sub_seurat_clusters == x]
    b.pv <- t.test(a, b, alternative = "two.sided", paired = TRUE)$p.value
    c.pv <- t.test(a, c, alternative = "two.sided", paired = TRUE)$p.value
    return(c("7d" = b.pv, 
             "2m" = c.pv))
  }) 
  sig.df <- data.frame(time_point = rownames(sig.df)[row(sig.df)],
                       sub_seurat_clusters = colnames(sig.df)[col(sig.df)],
                       p.value = c(sig.df))
  
  db_b <- db %>%
    dplyr::select(subject, sub_seurat_clusters, rel_diff_B) %>%
    dplyr::rename(rel_diff = rel_diff_B) %>%
    dplyr::mutate(time_point = "7d")
  
  db_c <- db %>%
    dplyr::select(subject, sub_seurat_clusters, rel_diff_C) %>%
    dplyr::rename(rel_diff = rel_diff_C) %>%
    dplyr::mutate(time_point = "2m")
  
  db_summ <- bind_rows(db_b, db_c) %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = subC_levels, ordered = T),
                  time_point = factor(time_point, levels = c("7d", "2m"), ordered = T)) %>%
    dplyr::group_by(sub_seurat_clusters) %>% 
    dplyr::mutate(height = max(rel_diff) + .3 * sd(rel_diff)) %>%
    dplyr::ungroup()
  
  sig.df <- merge(sig.df, db_summ, by = c("sub_seurat_clusters", "time_point"), all.x = T) %>%
    dplyr::select(sub_seurat_clusters, time_point, p.value, height) 
  sig.df <- sig.df[!duplicated(sig.df), ] %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = subC_levels, ordered = T),
                  time_point = factor(time_point, levels = c("7d", "2m"), ordered = T))
  
  
  db_summ <- db_summ %>%
    dplyr::group_by(sub_seurat_clusters, time_point) %>%
    dplyr::mutate(sub_seurat_clusters = factor(sub_seurat_clusters, levels = subC_levels, ordered = T),
                  outlier = ifelse(is_outlier(rel_diff), as.character(unique(subject)), as.character(NA))) %>%
    dplyr::ungroup()
  
  set.seed(12345)
  p <- ggplot(db_summ,
              aes(x = time_point, y = rel_diff, color = time_point, label = outlier)) +
    getBaseTheme() +
    theme(axis.text = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 24, face = "bold")) +
    xlab("Time point") + ylab("Relative change") +
    scale_x_discrete(limits =  time_points[-1],
                     labels = c("7d" = "B",
                                "2m" = "C")) +
    scale_y_continuous(expand = expansion(mult = c(0.04, .125))) +
    geom_hline(yintercept = 0, size = 0.4, linetype = 2) +
    geom_boxplot(width = 0.35, show.legend = F, outlier.size = 2, outlier.colour = NA) +
    geom_jitter(position = position_jitter(0.15), size = 2, show.legend = F) +
    geom_text_repel(size = 7, show.legend = F) +
    geom_text(data = sig.df,
              aes(x = time_point, y = height, label = ifelse(p.value < 0.05, "*", "")),
              fontface = "bold", size = 12, show.legend = F) +
    scale_fill_manual(name="", values = tp_colrs) +
    scale_color_manual(name="", values = tp_colrs) +
    facet_wrap(. ~ sub_seurat_clusters, 
               scales = "free", nrow = 1)
  plot(p)
  # ===================================
  id_levels <- as.vector(sapply(1:length(sub_lng), function(x){
    paste(sub_lng[x], time_points, sep=".")
  }))
  db <- db_plot %>%
    dplyr::mutate(id = paste(subject, time_point, sep=".")) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(frac = cell_cnt/sum(cell_cnt)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = factor(id, levels = id_levels, ordered = T))
  
  p <- ggplot() +
    getBaseTheme() +
    theme(axis.text.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 22, face = "bold"),
          strip.text = element_text(size = 22, face = "bold"),
          legend.text = element_text(size = 14)) +
    xlab("Time point") + ylab("Proportion of cells") + 
    geom_bar(data = db,
             aes(x = time_point, y = frac, fill = sub_seurat_clusters),
             width = 0.85, stat = "identity", position = position_stack(reverse = TRUE), show.legend = FALSE) +
    scale_fill_manual(name = "",
                      values = setNames(subC_colrs, levels(sobj))) +
    scale_x_discrete(limits =  time_points,
                     labels = c("0d" = "A",
                                "7d" = "B",
                                "2m" = "C")) +
    guides(fill = guide_legend(override.aes = list(size = 6), nrow = 1)) +
    facet_grid(. ~ subject)
  plot(p)
}
# ===================================
# find markers
if (TRUE) {
  top_genes <- all_markers_bc %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    dplyr::filter(avg_log2FC > lfc.th, 
                  p_val_adj < pv.th, 
                  pct.1 > 0.25, 
                  !grepl(pattern = "\\.", gene)) %>%
    dplyr::group_by(sub_seurat_clusters) %>%
    dplyr::top_n(n = 15, wt = avg_log2FC) %>%
    dplyr::arrange(avg_log2FC) %>%
    dplyr::ungroup()
  genes <- c(unique(top_genes$gene))
  
  mt_data <- FetchData(sobj, c("sub_seurat_clusters")) %>%
    tibble::rownames_to_column(var="barcode")
  
  genes_df <- FetchData(sobj, genes, slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    dplyr::mutate(cpm = expm1(log1pcpm)) %>% # back from log1pCPM to CPM (cell normalized)
    dplyr::group_by(gene, sub_seurat_clusters) %>%
    dplyr::summarise(mean_cpm = mean(cpm)) %>%
    dplyr::ungroup()   
}
# ===================================
if (TRUE) {
  genes_mtx <- genes_df %>%
    tidyr::spread(., key = gene, value = mean_cpm)
  row_names <- genes_mtx$sub_seurat_clusters
  genes_mtx$sub_seurat_clusters <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  # ===================================
  colnames_sorted <- unique(genes_df$sub_seurat_clusters)
  colnames_sorted <- colnames_sorted[mixedorder(colnames_sorted)]
  genes_mtx <- genes_mtx[, colnames_sorted]
  # ===================================
  cols.cor <- cor(genes_mtx, use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between samples (columns)
  rows.cor <- cor(t(genes_mtx), use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation between rows (genes)
  # ===================================
  n_clust <- n_subC
  my_row_hm <- hclust(as.dist(1 - rows.cor), method = "complete")
  my_row_hm <- cutree(tree = my_row_hm, k = n_clust)
  my_row_hm <- data.frame(Cluster = my_row_hm) %>%
    dplyr::mutate(Cluster = case_when(Cluster == 5 ~ subC_levels[5],
                                      Cluster == 4 ~ subC_levels[3],
                                      Cluster == 3 ~ subC_levels[2],
                                      Cluster == 2 ~ subC_levels[1],
                                      Cluster == 1 ~ subC_levels[4]))
  # ===================================
  my_colour_hm <- list(time_point = tp_colrs,
                       Cluster = setNames(subC_colrs, 
                                          subC_levels))
  # ===================================
  p <- pheatmap(genes_mtx,
                scale = "row", # This will scale the data to Z-scores (by gene) automatically
                show_colnames = T,
                show_rownames = T,
                labels_col = colnames(genes_mtx),
                labels_row = rownames(genes_mtx),
                color = colorpanel(1000,"blue","white","red"),
                clustering_method = "complete",
                clustering_distance_cols = "euclidean", # or use "correlation" or "euclidean"
                clustering_distance_rows = as.dist(1 - rows.cor), # or use "correlation"
                cluster_cols = TRUE,
                cluster_rows = TRUE,
                fontsize_col = 12,
                fontsize_row = 9,
                cutree_rows = n_clust,
                cutree_cols = n_clust,
                treeheight_row = 0,
                treeheight_col = 20.0,
                fontsize_number = 1,
                silent = F,
                border_color = NA,
                annotation_row = my_row_hm,
                annotation_colors = my_colour_hm,
                legend = TRUE,
                annotation_legend = FALSE
  )
}
# ===================================
if (TRUE) {
  deg.df <- deg_all_bc %>%
    dplyr::filter(abs(avg_log2FC) > lfc.th, 
                  p_val_adj < pv.th, 
                  !grepl(pattern = "\\.", gene)) %>%
    dplyr::mutate(reg_type = ifelse(avg_log2FC > 0, "up", "down")) %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name"))
  sde_genes <- unique(deg.df$gene)
  num_degs <- length(sde_genes)
}
# ===================================
if (TRUE) {
  deg.summ <- deg.df %>%
    dplyr::group_by(sub_seurat_clusters, time_point) %>%
    dplyr::summarise(n_degs = length(unique(gene)),
                     n_degs_up = sum(avg_log2FC > 0),
                     n_degs_down = sum(avg_log2FC <= 0)) %>%
    dplyr::ungroup() %>%
    tidyr::complete(., 
                    sub_seurat_clusters = subC_levels, 
                    time_point = unique(time_point), 
                    fill = list(n_degs = 0,
                                n_degs_up = 0,
                                n_degs_down = 0)) %>%
    tidyr::gather(., key = reg_type, value = count, -time_point, -sub_seurat_clusters) %>%
    dplyr::mutate(time_point = factor(time_point, levels = c("7d", "2m"), ordered = T),
                  sub_seurat_clusters = factor(sub_seurat_clusters, levels = subC_levels, ordered = T))
  
  db <- filter(deg.summ, reg_type != "n_degs") %>%
    dplyr::mutate(count = ifelse(reg_type == "n_degs_up", count, -count),
                  vjust = ifelse(count > 0, -0.5, 1.25))
  pd <- position_dodge2(width = 0.75, preserve = "single")
  p <- ggplot(db,
              aes(x = sub_seurat_clusters, y = count, fill = reg_type, label = abs(count), vjust = vjust)) +
    getBaseTheme() +
    theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          strip.text = element_text(size = 22, face = "bold"),
          legend.text = element_text(size = 14)
    ) +
    xlab("") + ylab("Number of sDEGs") +
    scale_y_continuous(labels = abs,
                       breaks = seq(-200, 200, 100)) +
    coord_cartesian(ylim = c(-200, 200)) +
    geom_bar(width = 0.75, stat = "identity", position = pd) +
    geom_text(data = filter(db, abs(count) > 0),
              position = pd, size = 5) +
    geom_hline(yintercept = 0, linetype = 1, size = 0.45) +
    scale_fill_manual(name = "",
                      labels = c("n_degs_up" = "Up-regulated sDEG", 
                                 "n_degs_down" = "Down-regulated sDEG"),
                      values = setNames(c("#e31a1c", "#1f78b4"), c("n_degs_up", "n_degs_down"))) +
    facet_grid(. ~ time_point, scales = "free",
               labeller = labeller(time_point = setNames(c("Time point B", "Time point C"),
                                                         c ("7d", "2m"))))
  plot(p)
}
# ===================================
if (TRUE) {
  p <- DimPlot(sobj,
                reduction = "umap",
                group.by = "sub_seurat_clusters",
                cols = setNames(subC_colrs, levels(sobj)),
                label = FALSE, 
                repel = TRUE,
                pt.size = 0.1,
                label.size = 6) + 
    getBaseTheme() +
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    guides(color = FALSE)
  plot(p)
}
# ===================================
if (TRUE) {
  genes <- c("CD69", # activation markers
             "IGHM", "IGHD", "IGHE", "IGHG1", 
             "IGHA1", "CD27", # memory BC  
             "CXCR4","CXCR5", "FCRL2", "FCRL5", 
             "IL10RA", 
             "CCR7", # B cells activated by antigen in the follicles upregulate CCR7
             "CCR9", "ISG15", "IFI44L", "MX1", "TCF4",
             "TMSB10", "TCL1A", "CD9", "CD19", "CD24", "CD38", 
    "HLA-A", "HLA-B", "HLA-C") 
  p <- DotPlot(sobj,
               features = genes, 
               cols = c("white", "red"),
               dot.min = 0.05,
               cluster.idents = TRUE) + 
    getBaseTheme() +
    RotatedAxis() +
    xlab("Gene") + ylab("SubClusters") +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 13, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"),
          legend.title=element_text(size=12, face = "bold"),
          legend.text=element_text(size=11, face = "bold"),
          legend.key.height=grid::unit(7.5, "points"), 
          legend.key.width=grid::unit(15, "points"),
          legend.spacing = grid::unit(40, "points")) +
    guides(fill = guide_legend(nrow = 1))
  plot(p)
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
    name = 'interferon'
  )
  # identifying highly correlated genes with module score
  scrs <- as.numeric(sobj_sc@meta.data$interferon1)
  correlations <- apply(as.matrix(sobj_sc@assays$RNA@data), 1, 
                        function(x){ cor(scrs, x) })
  top_ifn_genes_scrs <- names(head(sort(correlations[genes_ifn], decreasing = T), 10)) # run this twice for top 10 by 1:5 and 6:10
}
# ===================================
gene_overlap_test <- gene_overlap_fun(subC_levels, deg.df, mod)
# ===================================
if (TRUE) {
  db_plot <- gene_overlap_test %>%
    dplyr::filter(p_val < 0.05, signiture != "TBD") %>%
    dplyr::mutate(odds_ratio = ifelse(is.infinite(odds_ratio), max(odds_ratio[!is.infinite(odds_ratio)]), odds_ratio),
                  wrap_id = paste(time_point, reg_type),
                  wrap_id = factor(wrap_id, levels = c("7d down", "7d up", "2m down", "2m up"), ordered = T),
                  id = paste(signiture, module, sep = "_"),
                  id = factor(id, levels = unique(id[mixedorder(id)]), ordered = T))
  
  mds <- levels(db_plot$id)
  sec.axis.labels <- sapply(mds, function(x){
    unique(as.character(db_plot$module[db_plot$id == x]))
  }, USE.NAMES = FALSE)
  
  fst.axis.labels <- sapply(mds, function(x){
    unique(as.character(db_plot$signiture[db_plot$id == x]))
  }, USE.NAMES = FALSE)
  
  dotplot <- ggplot(db_plot,
                    aes(x = sub_seurat_clusters, y = id, fill = log10(odds_ratio), size = jaccard_index)) +   # size = `% Expressing`
    getBaseTheme() + 
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.title = element_text(size = 15, face = "bold"),
          legend.key.width = grid::unit(15, "points"),
          legend.key.height = grid::unit(20, "points"),
          legend.spacing = grid::unit(15, "points"),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(size = 15, face = "bold", angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 20, face = "bold")) +
    geom_point(shape = 21, color = "black") + 
    scale_size(name = "Jaccard\nIndex",
               range = c(1,7)) +
    xlab("Cluster") + ylab("Module Type") +
    scale_fill_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu")),
                         name = 'Odds Ratio\n[log10]') +
    scale_y_discrete(limit = mds,
                     labels = fst.axis.labels) +
    guides(y.sec = guide_axis_label_trans(label_trans = ~paste0(sec.axis.labels),
                                          "Module ID", angle = 0)) + 
    facet_wrap(. ~ wrap_id, nrow = 1,
               labeller = labeller(wrap_id = setNames(c("C-down", "C-up", "B-down", "B-up"),
                                                      c("2m down", "2m up", "7d down", "7d up"))))
  plot(dotplot)
}
# ===================================
if (TRUE) {
  genes <- unique(deg.df$gene[deg.df$gene %in% genes_ifn])
  mt_data <- FetchData(sobj, c("sub_seurat_clusters", "time_point")) %>%
    tibble::rownames_to_column(var="barcode")
  
  genes_df <- FetchData(sobj, genes, slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    dplyr::mutate(cpm = expm1(log1pcpm)) %>% 
    dplyr::group_by(sub_seurat_clusters, time_point, gene) %>%
    dplyr::summarise(mean_cpm = mean(cpm)) %>%
    dplyr::ungroup() 
  # ===================================
  genes_mtx <- genes_df %>%
    dplyr::mutate(id = paste(sub_seurat_clusters, time_point, sep = "_")) %>%
    dplyr::select(id, gene, mean_cpm) %>%
    tidyr::spread(., key = gene, value = mean_cpm) 
  row_names <- genes_mtx$id
  genes_mtx$id <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  # ===================================
  colnames_sorted <- unique(genes_df$sub_seurat_clusters)
  colnames_sorted <- colnames_sorted[mixedorder(colnames_sorted)]
  colnames_sorted <- paste(rep(colnames_sorted, each = length(time_points[1:3])), time_points[1:3], sep = "_")
  genes_mtx <- genes_mtx[rownames(genes_mtx)[mixedorder(rownames(genes_mtx))], colnames_sorted]
  # ===================================
  str_mtx <- matrix("", nrow = nrow(genes_mtx), ncol = ncol(genes_mtx), 
                    dimnames = list(rownames(genes_mtx), colnames(genes_mtx)))
  for(x in rownames(genes_mtx)) {
    db <- filter(deg.df, gene == x)
    str_mtx[x, paste(db$sub_seurat_clusters, db$time_point, sep = "_")] <- as.character(db$reg_type)
  }
  str_mtx[str_mtx == "up"] <- "u"
  str_mtx[str_mtx == "down"] <- "d"
  # ===================================
  my_col_hm <- data.frame(Cluster = unlist(lapply(strsplit(colnames_sorted, split = "_"), function(x) paste(x[1],x[2], sep = "_"))),
                          timePoint = unlist(lapply(strsplit(colnames_sorted, split = "_"), function(x) x[3]))) 
  rownames(my_col_hm) <- colnames_sorted
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
  # ===================================
  my_colour_hm <- list(timePoint = tp_colrs[1:3],
                       Cluster = setNames(subC_colrs, 
                                          unique(my_col_hm$Cluster)),
                       signiture = setNames(c("#fdc086", "#386cb0"),
                                            levels(my_row_hm$signiture)),
                       Module = setNames(brewer.pal(n = length(unique(my_row_hm$Module)), name = "Paired"),
                                         levels(my_row_hm$Module)),
                       Type = setNames(c("#fdc086", "#386cb0", "#f0027f", "grey70"),
                                       levels(my_row_hm$Type))
  )
  # ===================================
  bld_genes <- rownames(str_mtx)[rowSums(str_mtx[, grepl("7d", colnames(str_mtx))] == "u") == n_subC]
  p <- pheatmap(genes_mtx,
                scale = "row", # This will scale the data to Z-scores (by gene) automatically
                show_colnames = T,
                show_rownames = T,
                labels_col = my_col_hm$Cluster,
                labels_row = make_bold_names(genes_mtx, rownames, bld_genes),
                color = colorpanel(1000,"blue","white","red"),
                clustering_method = "complete",
                clustering_distance_cols = as.dist(1 - cols.cor), # or use "correlation"
                clustering_distance_rows = as.dist(1 - rows.cor), # or use "correlation"
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                fontsize_row = 10,
                fontsize_col = 10,
                gaps_col = seq(3, 12, 3),
                treeheight_row = 25,
                display_numbers = str_mtx,
                fontsize_number = 10,
                silent = F,
                border_color = "grey70",
                fontsize = 10,
                annotation_col = my_col_hm,
                annotation_row = my_row_hm[,c("Module", "Type"), drop = FALSE],
                annotation_colors = my_colour_hm
  )
}
# ===================================
library("viridis")
if (TRUE) {
  tmp <- data.frame(DC1 = sds@reducedDim[, "UMAP_1"],
                    DC2 = sds@reducedDim[, "UMAP_2"],
                    Timepoint = sobj@meta.data$sub_seurat_clusters[match(rownames(sds@clusterLabels), rownames(sobj@meta.data))])
  my_color <- setNames(subC_colrs, subC_levels)
  par(mfrow = c(1, 1), mai=c(1.0,1.0,0.25,0.25))
  plot(tmp$DC1, tmp$DC2,
       main = "",
       col = my_color[as.character(tmp$Timepoint)],
       xlab="", ylab="",
       font.lab = 2, font.main = 2,
       pch=16, cex = 0.4, axes = F)
  lines(sds@curves$curve1, 
        lwd = 5, col = 'black', cex = 0, lty = 1)
  lines(sds@curves$curve2, 
        lwd = 5, col = 'black', cex = 0, lty = 2)
  lines(sds@curves$curve3, 
        lwd = 5, col = 'black', cex = 0, lty = 6)
}
# ===================================
if (TRUE) {
  gg_titles <- c("Interferon")
  x_lim <- list(c(-1,2))
  gg_means <- list(gg_means <- sobj_sc@meta.data %>%
                     dplyr::group_by(sub_seurat_clusters) %>%
                     dplyr::summarise(means = mean(interferon1)))
  p <- RidgePlot(sobj_sc, 
                 features = c("interferon1"), 
                 assay = "RNA", slot = "data",
                 cols = setNames(subC_colrs, subC_levels),
                 sort = "increasing", same.y.lims = TRUE, combine = FALSE)
  p <- lapply(1:length(p), function(i) {
    p[[i]]  +  getBaseTheme() +
      theme(plot.title = element_blank(),
            strip.text = element_blank(),
            axis.title.x = element_text(size = 35, face = "bold"),
            axis.title.y = element_blank(),
            axis.text = element_text(size = 30, face = "bold"),
            legend.position='bottom',
            legend.title=element_blank(),
            legend.text=element_text(size=7),
            legend.key.height=grid::unit(10, "points"), 
            legend.key.width=grid::unit(20, "points"),
            legend.spacing = grid::unit(40, "points")) +
      ggtitle(gg_titles[i]) +
      xlab("Score") +
      coord_cartesian(xlim = c(x_lim[[i]][1], x_lim[[i]][2]), ylim = c(0,9)) +
      geom_vline(data = gg_means[[i]],
                 aes(xintercept = means, color = sub_seurat_clusters),
                 linetype = 1, size = 1.75) +
      scale_color_manual(name = "",
                         values = subC_colrs) +
      guides(fill = FALSE,
             color = FALSE)
  })
  p <- plot_grid(plotlist = p, nrow = 1)
  plot(p)
}
# ===================================
if (TRUE) {
  db_plot <- sobj_sc@meta.data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::select(barcode, sub_seurat_clusters, time_point, interferon1) %>%
    tidyr::gather(., key = module, value = score, -barcode, -sub_seurat_clusters, -time_point) %>%
    dplyr::mutate(time_point = factor(time_point, levels = time_points, ordered = T),
                  module = case_when(module == "interferon1" ~ "Interferon"),
                  module = factor(module, levels = unique(module), ordered = T)) %>%
    dplyr::group_by(module, sub_seurat_clusters, time_point) %>%
    dplyr::mutate(outlier = is_outlier(score)) %>%
    dplyr::ungroup()
  
  db_label <- db_plot %>%
    dplyr::group_by(module) %>%
    dplyr::mutate(y = max(score)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(module, sub_seurat_clusters, time_point) %>%
    dplyr::summarise(count = n(),
                     y = unique(y) + 0.1*unique(y)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time_point = factor(time_point, levels = time_points, ordered = T))
  error_bar <- db_plot %>%
    dplyr::group_by(module, sub_seurat_clusters, time_point) %>%
    dplyr::summarise(mean_sc = mean(score, na.rm = T),
                     sd_sc = sd(score, na.rm = T)) %>%
    dplyr::ungroup()
  
  pd <- position_dodge(width = 0.90)
  p <- ggplot(db_plot, 
              aes(x = sub_seurat_clusters, y = score, color = time_point, fill = time_point)) +
    getBaseTheme() +
    theme(plot.title = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 18, face = "bold"),
          strip.text = element_blank(),
          legend.text=element_text(size=12)) +
    xlab("") + ylab("Score") +
    geom_errorbar(data = error_bar, 
                  aes(x = sub_seurat_clusters, y = mean_sc,
                      ymin = mean_sc-sd_sc/2, ymax = mean_sc+sd_sc/2),
                  width=0.70, alpha=0.9, size=0.45, position = pd, show.legend = F) +
    stat_summary(fun = mean, geom = "crossbar", size = 0.6, width = 0.7, position = pd, show.legend = F) +
    geom_vline(xintercept = seq(1.5, 4.5, 1), linetype = 2, size= 0.35, alpha = 0.65) +
    scale_color_manual(name = "",
                       values = tp_colrs) +
    scale_fill_manual(name = "",
                      values = tp_colrs) +
    guides(color = guide_legend(override.aes = list(size = 0.75))) +
    facet_wrap(. ~ module, scale = "free")
  plot(p)
}
# ===================================
if (TRUE) {
  scrs <- c("interferon1")
  db_plot <- data.frame()
  df_smooth_all <- data.frame()
  for (scr in scrs) {
    mt_data <- FetchData(sobj_sc, c("sub_seurat_clusters", scr)) %>%
      tibble::rownames_to_column(var="barcode")  
    
    pt <- slingPseudotime(sds)
    db_c1 <- as.data.frame(pt) %>%
      tibble::rownames_to_column(var="barcode") %>%
      base::merge(., mt_data,
                  by = c("barcode"), all = TRUE) %>%
      dplyr::rename(pseudotime = curve1,
                    score = `scr`) %>%
      dplyr::mutate(trajectory = "Trajectory-1",
                    type = scr) %>%
      dplyr::select(barcode, sub_seurat_clusters, pseudotime, score, trajectory, type) %>%
      dplyr::mutate(pseudotime = scales::rescale(pseudotime, to=c(0,1)))
    model <- loess(score ~ pseudotime, data = db_c1)
    xseq <- seq(from=0, to=1, by=0.002)
    pred_1 <- predict(model, newdata = data.frame(pseudotime = xseq), se = FALSE)
    names(pred_1) <- xseq
    
    db_c2 <- as.data.frame(pt) %>%
      tibble::rownames_to_column(var="barcode") %>%
      base::merge(., mt_data,
                  by = c("barcode"), all = TRUE) %>%
      dplyr::rename(pseudotime = curve2,
                    score = `scr`) %>%
      dplyr::mutate(trajectory = "Trajectory-2",
                    type = scr) %>%
      dplyr::select(barcode, sub_seurat_clusters, pseudotime, score, trajectory, type) %>%
      dplyr::mutate(pseudotime = scales::rescale(pseudotime, to=c(0,1)))
    model <- loess(score ~ pseudotime, data = db_c2)
    xseq <- seq(from=0, to=1, by=0.002)
    pred_2 <- predict(model, newdata = data.frame(pseudotime = xseq), se = FALSE)
    names(pred_2) <- xseq
    
    db_c3 <- as.data.frame(pt) %>%
      tibble::rownames_to_column(var="barcode") %>%
      base::merge(., mt_data,
                  by = c("barcode"), all = TRUE) %>%
      dplyr::rename(pseudotime = curve3,
                    score = `scr`) %>%
      dplyr::mutate(trajectory = "Trajectory-3",
                    type = scr) %>%
      dplyr::select(barcode, sub_seurat_clusters, pseudotime, score, trajectory, type) %>%
      dplyr::mutate(pseudotime = scales::rescale(pseudotime, to=c(0,1)))
    model <- loess(score ~ pseudotime, data = db_c3)
    xseq <- seq(from=0, to=1, by=0.002)
    pred_3 <- predict(model, newdata = data.frame(pseudotime = xseq), se = FALSE)
    names(pred_3) <- xseq
    
    db_plot <- db_plot %>%
      dplyr::bind_rows(db_c1, db_c2, db_c3) %>%
      dplyr::filter(!is.na(pseudotime))
    
    df_smooth <- dplyr::bind_rows(data.frame(pseudotime = as.numeric(names(pred_1)),
                                             smmoth = as.vector(pred_1),
                                             trajectory = "Trajectory-1"),
                                  data.frame(pseudotime = as.numeric(names(pred_2)),
                                             smmoth = as.vector(pred_2),
                                             trajectory = "Trajectory-2"),
                                  data.frame(pseudotime = as.numeric(names(pred_3)),
                                             smmoth = as.vector(pred_3),
                                             trajectory = "Trajectory-3")) %>%
      dplyr::mutate(type = scr)
    
    df_smooth_all <- df_smooth_all %>%
      dplyr::bind_rows(df_smooth) 
  }
  db_plot <- db_plot %>%
    dplyr::mutate(type = factor(type, levels = scrs, ordered = T))
  
  df_smooth_all <- df_smooth_all %>%
    dplyr::mutate(smmoth = scales::rescale(smmoth, to=c(0,1)),
                  type = factor(type, levels = scrs, ordered = T))
  
  p <- ggplot() +
    getBaseTheme() +
    theme(strip.text = element_blank(),
          panel.spacing = unit(1, "lines"),
          axis.title.y = element_text(size = 19, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 10)),
          axis.title.y.right = element_text(size = 19, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 15)),
          axis.title.x = element_text(size = 19, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          legend.text = element_text(size=17),
          legend.key.height=grid::unit(10, "points"), 
          legend.key.width=grid::unit(40, "points"),
          legend.spacing = grid::unit(40, "points")) +
    xlab("Pseudotime (scaled)") + ylab("Score (scaled)") +
    geom_density(data = filter(db_plot, !(sub_seurat_clusters == "BC_2" & pseudotime < 0.1)),
                 aes(x = pseudotime, after_stat(count), fill = sub_seurat_clusters),
                 position = "fill", show.legend = F) +
    geom_line(data = df_smooth_all,
              aes(x = pseudotime, y = smmoth, linetype = trajectory),
              size = 1.0) +
    scale_y_continuous(sec.axis = sec_axis(trans=~., 
                                           name = "Cluster density (scaled)")) +
    scale_fill_manual(name = "",
                      values = setNames(subC_colrs,
                                        subC_levels)) +
    scale_linetype_manual(name = "",
                          values = setNames(c(1, 2, 6),
                                            c("Trajectory-1", "Trajectory-2", "Trajectory-3"))) +
    facet_grid(. ~ type, scales = "free",
               labeller = labeller(type = setNames(c("Interferon"), scrs)))
  plot(p)
}
# ===================================
if (TRUE) {
  mt_data <- FetchData(sobj, c("sub_seurat_clusters")) %>%
    tibble::rownames_to_column(var="barcode")
  
  pt <- as.data.frame(slingPseudotime(sds)) %>%
    tibble::rownames_to_column(var="barcode") %>%
    dplyr::rename(Trajectory_1 = curve1,
                  Trajectory_2 = curve2,
                  Trajectory_3 = curve3) %>%
    tidyr::gather(., key = trajectory, value = pseudotime, -barcode) 
  
  genes <- top_ifn_genes_scrs 
  
  db_plot <- FetchData(sobj, genes, slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    base::merge(., pt,
                by = c("barcode"), all = TRUE) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::group_by(gene, trajectory) %>%
    dplyr::mutate(gene = factor(gene, levels = genes, ordered = T), 
                  pseudotime = scales::rescale(pseudotime, to=c(0,1))) 
  
  df_smooth_all <- data.frame()
  for (gn in genes) {
    db <- filter(db_plot, trajectory == "Trajectory_1", gene == gn)
    model <- loess(log1pcpm ~ pseudotime, data = db)
    xseq <- seq(from=0, to=1, by=0.002)
    pred_1 <- predict(model, newdata = data.frame(pseudotime = xseq), se = FALSE)
    names(pred_1) <- xseq
    
    db <- filter(db_plot, trajectory == "Trajectory_2", gene == gn)
    model <- loess(log1pcpm ~ pseudotime, data = db)
    xseq <- seq(from=0, to=1, by=0.002)
    pred_2 <- predict(model, newdata = data.frame(pseudotime = xseq), se = FALSE)
    names(pred_2) <- xseq
    
    db <- filter(db_plot, trajectory == "Trajectory_3", gene == gn)
    model <- loess(log1pcpm ~ pseudotime, data = db)
    xseq <- seq(from=0, to=1, by=0.002)
    pred_3 <- predict(model, newdata = data.frame(pseudotime = xseq), se = FALSE)
    names(pred_3) <- xseq
    
    df_smooth <- dplyr::bind_rows(data.frame(pseudotime = as.numeric(names(pred_1)),
                                             smmoth = as.vector(pred_1),
                                             trajectory = "Trajectory-1"),
                                  data.frame(pseudotime = as.numeric(names(pred_2)),
                                             smmoth = as.vector(pred_2),
                                             trajectory = "Trajectory-2"),
                                  data.frame(pseudotime = as.numeric(names(pred_3)),
                                             smmoth = as.vector(pred_3),
                                             trajectory = "Trajectory-3")) %>%
      dplyr::mutate(gene = gn)
    
    df_smooth_all <- df_smooth_all %>%
      dplyr::bind_rows(df_smooth)
  }
  df_smooth_all <- df_smooth_all %>%
    dplyr::mutate(smmoth = scales::rescale(smmoth, to=c(0,1)),
                  gene = factor(gene, levels = genes, ordered = T))
  p <- ggplot() +
    getBaseTheme() +
    theme(strip.text = element_text(size = 24, face = "bold"),
          panel.spacing = unit(2, "lines"),
          axis.title.y = element_text(size = 19, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 10)),
          axis.title.y.right = element_text(size = 19, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 15)),
          axis.title.x = element_text(size = 19, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          legend.text = element_text(size=17),
          legend.key.height=grid::unit(10, "points"), 
          legend.key.width=grid::unit(30, "points"),
          legend.spacing = grid::unit(40, "points")) +
    xlab("Pseudotime (scaled)") + ylab("Expression level (scaled)") +
    geom_density(data = filter(db_plot, !(sub_seurat_clusters == "BC_2" & pseudotime < 0.1), gene %in% top_ifn_genes_scrs[6:10]),
                 aes(x = pseudotime, after_stat(count), fill = sub_seurat_clusters),
                 position = "fill", show.legend = F) +
    geom_line(data = filter(df_smooth_all, gene %in% top_ifn_genes_scrs[6:10]),
              aes(x = pseudotime, y = smmoth, linetype = trajectory),
              size = 1.0) +
    scale_y_continuous(sec.axis = sec_axis(trans=~., 
                                           name = "Cluster density (scaled)")) +
    scale_fill_manual(name = "",
                      values = setNames(subC_colrs,
                                        subC_levels)) +
    scale_linetype_manual(name = "",
                          values = setNames(c(1, 2, 6),
                                            c("Trajectory-1", "Trajectory-2", "Trajectory-3"))) +
    facet_grid(. ~ gene, scales = "free")
  plot(p)
}
# ===================================
if (TRUE) {
  genes_ls <- list()
  for (subC in subC_levels) {
    for (tp in time_points[2:3]) {
      for (regtype in c("up", "down")) {
        genes <- deg.df %>%
          dplyr::filter(reg_type == regtype, sub_seurat_clusters == subC, time_point == tp) %>%
          dplyr::pull(gene) 
        genes_ls[[length(genes_ls) + 1]] <- unique(genes)
        names(genes_ls)[length(genes_ls)] <- paste(subC, tp, regtype, sep = "_")
      }
    }
  }
  
  for (regtype in c("up", "down")) {
    for (tp in time_points[2:3]) {
      genes_subs_ls <- genes_ls[grep(pattern = paste(tp, regtype, sep = "_"), names(genes_ls))]
      x <- unlist(genes_subs_ls, use.names = F)
      print(x[grep(pattern = "^MT-", x)])
      p <- ggVennDiagram(genes_subs_ls, 
                         category.names = substr(names(genes_subs_ls), 1, 4),
                         set_size = 0,
                         set_color = setNames(subC_colrs,
                                              subC_levels),
                         label = "count",
                         label_size = 12,
                         label_alpha = 0) +
        scale_fill_gradientn(colours = brewer.pal(n = 9, name = ifelse(regtype == "up", "Reds", "Blues"))[-c(7:9)]) +
        scale_color_manual(values = setNames(rep("black", n_subC),
                                             subC_levels)) +
        guides(fill = FALSE)
      plot(p)
    }  
  }
}
# ===================================
if (TRUE) {
  mt_data <- FetchData(sobj, c("sub_seurat_clusters", "time_point")) %>%
    tibble::rownames_to_column(var="barcode")
  
  pt <- slingPseudotime(sds)
  db_plot <- as.data.frame(pt) %>%
    tibble::rownames_to_column(var="barcode") %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    tidyr::gather(., key = "lineage", value = "pseudotime", -barcode, -sub_seurat_clusters, -time_point) %>%
    dplyr::mutate(time_point = factor(time_point, levels = time_points, ordered = T)) %>%
    dplyr::group_by(lineage) %>%
    dplyr::mutate(pseudotime = scales::rescale(pseudotime, to=c(0,1))) %>%
    dplyr::ungroup()

  pd <- position_dodge(width = 0.90)
  p <- ggplot(data = filter(db_plot, pseudotime > 0),
              aes(x = pseudotime, y = time_point, fill = sub_seurat_clusters, color = sub_seurat_clusters)) +
    getBaseTheme() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text=element_text(size=12)) +
    xlab("Pseudotime (scaled)") + ylab("Time point") +
    geom_boxplot(position = pd, width = 0.60, outlier.size = 0.6, show.legend = F, outlier.colour = NA, na.rm = T) +
    scale_y_discrete(limits =  c(time_points),
                     labels = c("0d" = "A",
                                "7d" = "B",
                                "2m" = "C")) +
    scale_fill_manual(name = "",
                      values = setNames(subC_colrs,
                                        subC_levels)) +
    scale_color_manual(name = "",
                       values = setNames(subC_colrs,
                                         subC_levels)) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 4)),
           fill= guide_legend(nrow = 1, override.aes = list(size = 4)))
  plot(p)
}
# ===================================