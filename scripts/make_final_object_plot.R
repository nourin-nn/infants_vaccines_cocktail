library("dplyr")
library("Seurat")
library("pheatmap")
library("ggplot2")
library("gplots")
library("grid")
library("gridExtra")
library("ggbeeswarm")
library("cowplot")
library("ggrepel")
library("tidyverse")
library("ggdendro")
library("ggtree")
library("patchwork") 
library("aplot")
library("treemapify")
# ===================================
sobj_res <- readRDS(file.path(analysis_path, "final_object.rds"))
sobj <- sobj_res
# ===================================
Idents(sobj) <- "cell_type"
# ===================================
reads_per_cell <- colSums(sobj@assays$RNA@counts)
mean(reads_per_cell) # reads/cell
sd(reads_per_cell)

genes_per_cell <- colSums(sobj@assays$RNA@counts != 0)
mean(genes_per_cell) # genes/cell
sd(genes_per_cell)
# ===================================
if (TRUE) {
  cell_type_annot_vec <- unique(sobj@meta.data$cell_type)  
  n_cell_type_annot <- length(cell_type_annot_vec)
  
  cell_type_color <- c()
  
  sort(cell_type_annot_vec[grepl("CD4", cell_type_annot_vec)])
  cell_type_color <- c(cell_type_color, 
                       setNames(brewer.pal(n = 9, name = "Greens")[-c(1,2,3,9)], 
                                c("CD4 TC/C0", "CD4 TC/C1", "CD4 TC/C5", "CD4 TC/C7", "CD4 TC/C8"))) 
  
  sort(cell_type_annot_vec[grepl("CD8", cell_type_annot_vec)])
  cell_type_color <- c(cell_type_color, 
                       setNames(brewer.pal(n = 9, name = "Purples")[c(4,9,7)], 
                                c("CD8 TC/C2", "CD8 TC/C11", "CD8 TC/C10"))) 
  
  sort(cell_type_annot_vec[grepl("NK", cell_type_annot_vec)])
  cell_type_color <- c(cell_type_color, 
                       setNames(brewer.pal(n = 9, name = "Reds")[c(6,8)], 
                                c("NKC/C4", "NKC/C13"))) 
  
  sort(cell_type_annot_vec[grepl("PC|BC", cell_type_annot_vec)])
  cell_type_color <- c(cell_type_color, 
                       setNames(c("#df65b0", "#980043", "#c994c7"), 
                                c("BC/C3", "BC/C12", "PC/C17"))) 
  
  sort(cell_type_annot_vec[grepl("Mono|cDC", cell_type_annot_vec)])
  cell_type_color <- c(cell_type_color, 
                       setNames(c("blue", "#9ecae1", "#3182bd"), 
                                c("CD16 Mono/C9", "CD14 Mono/C6","cDC/C16")))
  
  sort(cell_type_annot_vec[grepl("pDC", cell_type_annot_vec)])
  cell_type_color <- c(cell_type_color, 
                       setNames(c("black"), 
                                c("pDC/C15"))) 
  
  sort(cell_type_annot_vec[grepl("Mgk|Eryt", cell_type_annot_vec)])
  cell_type_color <- c(cell_type_color, 
                       setNames(c("#e7298a",alpha("#980043", 0.75)), 
                                c("Eryt/C14", "Mgk/C18"))) 
  
  col_order <- gsub(".*/", "", names(cell_type_color)) 
}
# ===================================
if (TRUE) {
  genes_ls <- list(
    "CD4 TC" = c("CD3D", "IL7R"), 
    "CD8 TC" = c("CD8A", "CD8B"), 
    "NKC" = c("NKG7", "GNLY"), 
    "BC" = c("CD79A", "MS4A1"), 
    "PC" = c("IGHA1", "IGHG1", "JCHAIN"), 
    "Mono" = c("MS4A7", "VCAN", "CD14", "FCN1"), 
    "cDC" = c("CD1C" ,"CST3", "FCER1A", "CLEC10A"), 
    "pDC" = c("TCF4", "IL3RA" ,"IRF7" ,"LILRA4"), 
    "Eryt" = c("HBA1" ,"HBA2", "HBB"), 
    "Mgk" = c("PF4", "PPBP")
    )
  # ===================================
  mt_data <- FetchData(sobj, c("cell_type")) %>%
    tibble::rownames_to_column(var="barcode")
  
  genes_df <- FetchData(sobj, unlist(genes_ls, use.names = F), slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    dplyr::mutate(cpm = expm1(log1pcpm), # back from log1pCPM to CPM (cell normalized)
                  cell_type_mrg = gsub("\\/.*", "", cell_type),
                  cell_type_mrg = ifelse(grepl(pattern = "Mono", cell_type_mrg), "Mono", cell_type_mrg)) %>% 
    dplyr::group_by(gene, cell_type_mrg) %>%
    dplyr::summarise(mean_cpm = mean(cpm)) %>%
    dplyr::ungroup()   
  # ===================================
  genes_mtx <- genes_df %>%
    dplyr::select(cell_type_mrg, gene, mean_cpm) %>%
    tidyr::spread(., key = gene, value = mean_cpm) 
  row_names <- genes_mtx$cell_type_mrg
  genes_mtx$cell_type_mrg <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  genes_mtx <- genes_mtx[rowSums(genes_mtx) != 0, ]
  genes_mtx <- genes_mtx[unlist(genes_ls, use.names = F), names(genes_ls)]
  # ===================================
  my_row_hm <- stack(x = genes_ls) %>%
    tibble::column_to_rownames(var="values") %>%
    dplyr::rename(Cluster = ind) 
  # ===================================
  my_colour_hm <- list(Cluster = setNames(c(cell_type_color["CD4 TC/C0"],cell_type_color["CD8 TC/C11"],cell_type_color["NKC/C4"],
                                            cell_type_color["BC/C3"],cell_type_color["PC/C17"],cell_type_color["CD16 Mono/C9"],
                                            cell_type_color["pDC/C15"],cell_type_color["cDC/C16"],cell_type_color["Eryt/C14"],
                                            cell_type_color["Mgk/C18"]), 
                                          names(genes_ls)))
  # ===================================
  p <- pheatmap(genes_mtx,
                scale = "row", # This will scale the data to Z-scores (by gene) automatically
                show_colnames = T,
                show_rownames = T,
                color = colorpanel(1000,"#6a3d9a","white","#941751"),
                clustering_method = "complete",
                clustering_distance_cols = "euclidean", # or use "correlation"
                clustering_distance_rows = "correlation", # or use "correlation"
                cluster_cols = F,
                cluster_rows = F,
                fontsize_col = 12,
                fontsize_row = 12,
                gaps_row = cumsum(lengths(genes_ls)),
                silent = F,
                border_color = "white",
                legend = TRUE,
                annotation_legend = FALSE,
                annotation_row = my_row_hm,
                annotation_colors = my_colour_hm
  )
}
# ===================================
if (TRUE) {
  cell_types <- unique(sobj@meta.data$cell_type)
  cltyp <- "CD4|CD8|NKC"
  clusts <- cell_types[grepl(pattern = cltyp, cell_types)]
  sobj_flt <- subset(sobj, subset = cell_type %in% clusts)

  genes <- c("CD3D", "CD4", "CD8A", "KLRF1", "NCAM1") 
  
  mt_data <- FetchData(sobj_flt, c("sub_seurat_clusters")) %>%
    tibble::rownames_to_column(var="barcode")
  
  genes_df <- FetchData(sobj_flt, genes, slot = "data") %>%
    tibble::rownames_to_column(var = "barcode") %>%
    tidyr::gather(., key = gene, value = log1pcpm, -barcode) %>%
    base::merge(., mt_data,
                by = c("barcode"), all = TRUE) %>%
    dplyr::mutate(cpm = expm1(log1pcpm)) %>% # back from log1pCPM to CPM (cell normalized)
    dplyr::group_by(gene, sub_seurat_clusters) %>%
    dplyr::summarise(mean_cpm = mean(cpm)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene = factor(gene, levels = genes, ordered = T))
  
  genes_mtx <- genes_df %>%
    dplyr::select(sub_seurat_clusters, gene, mean_cpm) %>%
    tidyr::spread(., key = gene, value = mean_cpm) 
  row_names <- genes_mtx$sub_seurat_clusters
  genes_mtx$sub_seurat_clusters <- NULL
  genes_mtx <- as.matrix(genes_mtx)
  rownames(genes_mtx) <- row_names
  genes_mtx <- t(genes_mtx)
  p <- pheatmap(genes_mtx,
                scale = "row", # This will scale the data to Z-scores (by gene) automatically
                show_colnames = T,
                show_rownames = T,
                color = colorpanel(1000,"#6a3d9a","white","#b2182b"),
                clustering_method = "complete",
                clustering_distance_cols = "euclidean", # or use "correlation"
                clustering_distance_rows = "correlation", # or use "correlation"
                cluster_cols = F,
                cluster_rows = F,
                fontsize_col = 12,
                fontsize_row = 12,
                silent = F,
                border_color = "white",
                legend = TRUE)
}
# ===================================
if (TRUE) {
  p <- DimPlot(sobj,
               reduction = "umap",
               group.by = "cell_type",
               cols = cell_type_color,
               label = FALSE, 
               repel = TRUE,
               label.size = 4) + 
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
  db_plot <- sobj@meta.data %>%
    dplyr::mutate(lbl = case_when(cell_type == "CD4 TC/C0" ~ "CD4_0",
                                  cell_type == "CD4 TC/C1" ~ "CD4_1",
                                  cell_type == "CD4 TC/C5" ~ "CD4_2",
                                  cell_type == "CD4 TC/C7" ~ "CD4_3",
                                  cell_type == "CD4 TC/C8" ~ "CD4_4",
                                  # ===================================
                                  cell_type == "CD8 TC/C2" ~ "CD8_0",
                                  cell_type == "CD8 TC/C10" ~ "CD8_1",
                                  cell_type == "CD8 TC/C11" ~ "CD8_2",
                                  # ===================================
                                  cell_type == "NKC/C4" ~ "NKC_0",
                                  cell_type == "NKC/C13" ~ "NKC_1",
                                  # ===================================
                                  cell_type == "BC/C3" ~ "BC_0",
                                  cell_type == "BC/C12" ~ "BC_1",
                                  # ===================================
                                  cell_type == "PC/C17" ~ "PC_0",
                                  # ===================================
                                  cell_type == "CD14 Mono/C6" ~ "Mono_0",
                                  cell_type == "CD16 Mono/C9" ~ "Mono_1",
                                  # ===================================
                                  cell_type == "cDC/C16" ~ "cDC_0",
                                  # ===================================
                                  cell_type == "pDC/C15" ~ "pDC_0",
                                  # ===================================
                                  cell_type == "Eryt/C14" ~ "Eryt_0",
                                  # ===================================
                                  cell_type == "Mgk/C18" ~ "Mgk_0")) 
  # ===================================
  lbls_level <- c("CD4_0", "CD4_1", "CD4_2", "CD4_3", "CD4_4",
                      "CD8_0", "CD8_1", "CD8_2",
                      "NKC_0", "NKC_1",
                      "BC_0", "BC_1",
                      "PC_0", 
                      "Mono_0", "Mono_1",
                      "cDC_0", "pDC_0", "Eryt_0", "Mgk_0", 
                      "All")
  # ===================================
  tbl <- table(sobj@meta.data$time_point)
  db_p <- db_plot %>%
    dplyr::group_by(lbl, time_point) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(lbl) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(data.frame(lbl = "All",
                                time_point = names(tbl),
                                count = as.vector(tbl),
                                frac = as.vector(tbl/sum(tbl)))) %>%
    dplyr::mutate(lbl = factor(lbl, levels = rev(lbls_level), ordered = T),
                  time_point = factor(time_point, levels = time_points, ordered = T))
  p <- ggplot() +
    getBaseTheme() +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15)) +
    xlab("") + ylab("Propotion of cells") + 
    geom_bar(data = db_p,
             aes(x = lbl, y = frac, fill = time_point),
             width = 0.85, stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(name = "",
                      labels = c("0d" = "A",
                                 "7d" = "B",
                                 "2m" = "C"),
                      values = tp_colrs) +
    coord_flip(ylim = c(0,1)) + 
    guides(fill = guide_legend(override.aes = list(size = 7.5)))
  plot(p)
  # ===================================
  tbl <- table(sobj@meta.data$batch)
  db_p <- db_plot %>%
    dplyr::group_by(lbl, batch) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(lbl) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(data.frame(lbl = "All",
                                batch = names(tbl),
                                count = as.vector(tbl),
                                frac = as.vector(tbl/sum(tbl)))) %>%
    dplyr::mutate(lbl = factor(lbl, levels = lbls_level, ordered = T),
                  batch = factor(batch, levels = unique(batch), ordered = T))
  p1 <- ggplot() +
    getBaseTheme() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 24, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 17, face = "bold", angle = 45, hjust = 1, vjust = 1),
          legend.position = "right",
          legend.title = element_text(size=20, face = "bold"),
          legend.text = element_text(size=20)) +
    xlab("") + ylab("Propotion of cells") + 
    geom_bar(data = db_p,
             aes(x = lbl, y = frac, fill = batch),
             width = 0.85, stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(name = "Batch",
                      values = batch_colrs) +
    guides(fill = guide_legend(override.aes = list(size = 7.5)))
  # ===================================
  tbl <- table(sobj@meta.data$sample)
  db_p <- db_plot %>%
    dplyr::group_by(lbl, sample) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(lbl) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(data.frame(lbl = "All",
                                sample = names(tbl),
                                count = as.vector(tbl),
                                frac = as.vector(tbl/sum(tbl)))) %>%
    dplyr::mutate(lbl = factor(lbl, levels = lbls_level, ordered = T),
                  sample_tmp = paste0("S", 1:length(unique(sample)))[as.factor(sample)],
                  sample_tmp = factor(sample_tmp, levels = unique(sample_tmp), ordered = T))
  
  sample_lng_colrs_tmp <- sample_lng_colrs
  names(sample_lng_colrs_tmp) <- unique(db_p$sample_tmp)
  
  p2 <- ggplot() +
    getBaseTheme() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 24, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 17, face = "bold", angle = 45, hjust = 1, vjust = 1),
          legend.position = "right",
          legend.title = element_text(size=20, face = "bold"),
          legend.text = element_text(size=20)) +
    xlab("") + ylab("Propotion of cells") + 
    geom_bar(data = db_p,
             aes(x = lbl, y = frac, fill = sample_tmp),
             width = 0.85, stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(name = "Sample",
                      values = sample_lng_colrs_tmp) +
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 7.5)))
  # ===================================
  tbl <- table(sobj@meta.data$subject)
  db_p <- db_plot %>%
    dplyr::group_by(lbl, subject) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(lbl) %>%
    dplyr::mutate(frac = count / sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(data.frame(lbl = "All",
                                subject = names(tbl),
                                count = as.vector(tbl),
                                frac = as.vector(tbl/sum(tbl)))) %>%
    dplyr::mutate(lbl = factor(lbl, levels = lbls_level, ordered = T),
                  subject = factor(subject, levels = unique(subject), ordered = T))
  p3 <- ggplot() +
    getBaseTheme() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 24, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 17, face = "bold", angle = 45, hjust = 1, vjust = 1),
          legend.position = "right",
          legend.title = element_text(size=20, face = "bold"),
          legend.text = element_text(size=20)) +
    xlab("") + ylab("Propotion of cells") + 
    geom_bar(data = db_p,
             aes(x = lbl, y = frac, fill = subject),
             width = 0.85, stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(name = "Infant",
                      values = sub_lng_colrs) +
    guides(fill = guide_legend(override.aes = list(size = 7.5)))
  # ===================================
  p <- plot_grid(p1,p2,p3, 
                 rel_widths = c(0.32,0.36,0.32),
                 nrow = 1)
  plot(p)
}
# ===================================
if (TRUE) {
  db_plot <- sobj@meta.data %>%
    dplyr::mutate(cell_type_mrg = gsub("\\/.*", "", cell_type),
                  cell_type_mrg = ifelse(grepl("Mono", cell_type_mrg), "Mono", cell_type_mrg)) %>%
    dplyr::group_by(cell_type_mrg) %>%
    dplyr::summarise(cell_type_cnt = n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cell_type_cnt) %>%
    dplyr::mutate(cell_type_mrg = factor(cell_type_mrg, levels = unique(cell_type_mrg), ordered = T))
  
  p <- ggplot() +
    getBaseTheme() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 19, face = "bold", angle = 30, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 20, face = "bold")) +
    xlab("") + ylab("Number of cells") + 
    geom_bar(data = db_plot,
             aes(x = cell_type_mrg, y = cell_type_cnt),
             width = 0.75, stat = "identity", fill = "steelblue") +
    geom_text(data = db_plot,
              aes(x = cell_type_mrg, y = cell_type_cnt, label = scales::comma(cell_type_cnt, accuracy = 1)),
              size = 5, nudge_y = 0.25) +
    scale_y_continuous(trans = scales::log10_trans(),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))
  plot(p)
}
# ===================================