## =========================================================
##  Pipeline DESeq2 — TEST (Anotação GTF)
##  Autor: Walter + ChatGPT
##  Saída principal de figuras: SVG + PDF
## =========================================================

suppressPackageStartupMessages({
  libs <- c("DESeq2","ggplot2","dplyr","readr","pheatmap",
            "RColorBrewer","matrixStats","stringr","svglite")
  for (p in libs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Pacote ausente: ", p, ". Instale e rode de novo.")
    }
  }
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
  library(stringr)
  library(svglite)
  
  has_xlsx <- requireNamespace("openxlsx", quietly = TRUE)
  if (has_xlsx) library(openxlsx)
})

## ============================
## 0) ENTRADAS & SAÍDAS
## ============================

## >>> AJUSTE AQUI (TEST) <<<
matriz_file <- "/Volumes/Expansion/doc_lastrun/test_tsv/matriz_kallisto.isoform.counts.matrix"
gtf_file    <- "/Users/walterfranco/Downloads/gtfmeso/ncbi_dataset/data/GCF_017639785.1/genomic.gtf"
out_root    <- "/Users/walterfranco/Library/CloudStorage/OneDrive-Pessoal(2)/Doutorado iec/Experimentos/Analise_seq/analisados/test/test_22_12"

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

## subpastas
mk <- function(...) {
  d <- file.path(...)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}
dir_logs  <- mk(out_root, "logs")
dir_tabs  <- mk(out_root, "tabelas")
dir_vol   <- mk(out_root, "volcano")
dir_pca   <- mk(out_root, "pca")
dir_pca_g <- mk(dir_pca, "global")
dir_pca_c <- mk(dir_pca, "por_condicao")
dir_hm    <- mk(out_root, "heatmap")
dir_sum   <- mk(out_root, "resumos")

## log
logfile <- file.path(dir_logs, "pipeline_DESeq2_TEST.log")
zz <- file(logfile, open = "wt")
sink(zz, type = "message")
message("## Log iniciado: ", as.character(Sys.time()))

## thresholds globais
padj_thr <- 0.05
lfc_thr  <- 1

## =========================================================
## (PARÂMETROS) DEGs + variância (para PCA/Heatmap)
## =========================================================
PCA_MAX_GENES_VAR       <- 3000  # PCA global/contraste: no máximo 3000 genes (se tiver muito)
HEATMAP_MAX_GENES_VAR   <- 800   # Heatmap: no máximo 800 genes (legibilidade)

## =========================================================
## 1) Funções auxiliares
## =========================================================

## 1.1 Leitura simples da matriz
read_counts_matrix <- function(path){
  message(">> Lendo matriz: ", path)
  if (!file.exists(path)) stop("Arquivo não encontrado: ", path)
  d <- read.table(path, header = TRUE, row.names = 1, sep = "\t",
                  check.names = FALSE, quote = "", comment.char = "")
  stopifnot(is.data.frame(d))
  message("   Encontrados: ", nrow(d), " features × ", ncol(d), " amostras")
  d
}

## 1.2 Inferir condição a partir do nome da amostra
## CN* -> CN ; 03DPI* -> 03dpi ; etc.
get_condition <- function(s) {
  s2 <- toupper(s)
  if (grepl("^CN", s2)) return("CN")
  m <- regmatches(s2, regexpr("^\\d{2}DPI", s2))
  if (length(m) == 1 && nchar(m) > 0) {
    return(gsub("DPI", "dpi", tolower(m)))
  }
  return("CN")
}

## 1.3 Limpar IDs dos transcritos para bater com GTF
clean_tx_id <- function(x){
  x <- sub("^rna-", "", x, ignore.case = TRUE)
  x <- sub("(_\\d+)$", "", x)
  x
}

## 1.4 Tipo do ID
tag_tipo <- function(ids) {
  ifelse(grepl("^rna-?XM_", ids, ignore.case = TRUE), "XM",
         ifelse(grepl("^MSTRG", ids, ignore.case = TRUE), "MSTRG", "Other"))
}

## 1.5 Shrink de log2FC (robusto; se falhar usa bruto)
do_shrink <- function(dds, contrast_vec) {
  if (requireNamespace("apeglm", quietly = TRUE)) {
    s <- try(lfcShrink(dds, contrast = contrast_vec, type = "apeglm"),
             silent = TRUE)
    if (!inherits(s, "try-error")) return(s)
  }
  if (requireNamespace("ashr", quietly = TRUE)) {
    s <- try(lfcShrink(dds, contrast = contrast_vec, type = "ashr"),
             silent = TRUE)
    if (!inherits(s, "try-error")) return(s)
  }
  message("(!) apeglm/ashr falharam ou ausentes — usando log2FoldChange bruto.")
  res <- results(dds, contrast = contrast_vec)
  res$log2FoldChange
}

## 1.6 Volcano plot (NatureStyle) — SVG + PDF
plot_volcano <- function(df, titulo, outfile_base,
                         padj_thr = 0.05, lfc_thr = 1) {
  df <- df %>%
    mutate(
      Signif = ifelse(!is.na(padj) &
                        padj <= padj_thr &
                        abs(log2FC_use) >= lfc_thr,
                      "Significativo","NS")
    )
  p <- ggplot(df, aes(x = log2FC_use, y = -log10(pvalue), color = Signif)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    scale_color_manual(values = c("NS" = "grey70", "Significativo" = "#D95F02")) +
    theme_minimal(base_size = 12) +
    labs(x = "log2 Fold Change (shrink)",
         y = "-log10(p-value)",
         title = titulo)
  
  ggsave(paste0(outfile_base, ".svg"), p, width = 7, height = 6, device = "svg")
  ggsave(paste0(outfile_base, ".pdf"), p, width = 7, height = 6)
}

## 1.7 PCA por contraste (condição vs CN) usando DEGs do contraste + var — SVG + PDF
plot_pca_contraste <- function(ct, data, coldata, outdir, degs_ct = NULL) {
  message("   - PCA: ", ct, " vs CN ...")
  keep <- coldata$condition %in% c("CN", ct)
  data_sub    <- data[, keep, drop = FALSE]
  coldata_sub <- droplevels(coldata[keep, , drop = FALSE])
  
  dds_sub <- DESeqDataSetFromMatrix(countData = data_sub,
                                    colData   = coldata_sub,
                                    design    = ~ condition)
  vsd_sub <- vst(dds_sub, blind = TRUE)
  
  if (!is.null(degs_ct)) {
    degs_ct <- intersect(degs_ct, rownames(vsd_sub))
    if (length(degs_ct) >= 2) {
      mat <- assay(vsd_sub)[degs_ct, , drop = FALSE]
      vv  <- matrixStats::rowVars(mat)
      if (length(vv) > PCA_MAX_GENES_VAR) {
        ord <- order(vv, decreasing = TRUE)
        degs_ct <- degs_ct[ord][seq_len(PCA_MAX_GENES_VAR)]
      }
      vsd_sub <- vsd_sub[degs_ct, ]
    }
  }
  
  pca_sub <- plotPCA(vsd_sub, intgroup = "condition",
                     ntop = min(nrow(vsd_sub), nrow(vsd_sub)),
                     returnData = TRUE)
  pv <- round(100 * attr(pca_sub, "percentVar"))
  
  p <- ggplot(pca_sub, aes(PC1, PC2, color = condition)) +
    geom_point(size = 10) +
    xlab(paste0("PC1: ", pv[1], "% var")) +
    ylab(paste0("PC2: ", pv[2], "% var")) +
    theme_minimal(base_size = 14) +
    ggtitle(paste0("PCA - ", ct, " vs CN (DEGs + var)"))
  
  ggsave(file.path(outdir, paste0("PCA_", ct, "_vs_CN.svg")),
         p, width = 7, height = 6, device = "svg")
  ggsave(file.path(outdir, paste0("PCA_", ct, "_vs_CN.pdf")),
         p, width = 7, height = 6)
  
  message("PCA salvo: ", ct, " vs CN")
}

## 1.8 Ler GTF e montar tabela (transcript_id -> gene_id, gene_name, biotype)
parse_gtf_tx2gene <- function(gtf_file){
  message(">> Lendo GTF para anotação: ", gtf_file)
  if (!file.exists(gtf_file))
    stop("GTF não encontrado: ", gtf_file)
  
  gtf <- read_tsv(
    gtf_file,
    comment = "#",
    col_names = c("seqname","source","feature","start","end",
                  "score","strand","frame","attribute"),
    col_types = "ccccccccc"
  )
  
  gtf_tx <- gtf %>% filter(feature %in% c("transcript","mRNA"))
  
  extract_attr <- function(x, key) {
    m <- str_match(x, paste0(key, " \"([^\"]+)\""))
    m[,2]
  }
  
  tx2gene <- gtf_tx %>%
    mutate(
      transcript_id = extract_attr(attribute, "transcript_id"),
      gene_id       = extract_attr(attribute, "gene_id"),
      gene_name     = extract_attr(attribute, "gene"),
      gene_biotype  = extract_attr(attribute, "gene_biotype")
    ) %>%
    select(transcript_id, gene_id, gene_name, gene_biotype) %>%
    distinct()
  
  message("   Transcritos anotados únicos: ", nrow(tx2gene))
  tx2gene
}

## =========================================================
## 2) Ler matriz & montar DESeq2
## =========================================================

data <- read_counts_matrix(matriz_file)

## garantir inteiros e filtro leve
data <- round(as.matrix(data))
data <- data[rowSums(data) >= 2, , drop = FALSE]

samples   <- colnames(data)
condition <- vapply(samples, get_condition, character(1))

## >>> AJUSTE: níveis típicos do TEST (03/05/10/15/30/CN)
condition <- factor(condition,
                    levels = c("03dpi","05dpi","10dpi","15dpi","30dpi","CN"))
coldata <- data.frame(condition = condition, row.names = samples)

stopifnot(ncol(data) == nrow(coldata))
message("Resumo coldata (n amostras por condição):")
print(table(coldata$condition))

message(">> Rodando DESeq2 (global) ...")
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData   = coldata,
                              design    = ~ condition)
dds <- DESeq(dds)

message(">> Transformação VST ...")
vsd <- vst(dds, blind = TRUE)

tx2gene <- parse_gtf_tx2gene(gtf_file)

if (has_xlsx) wb <- openxlsx::createWorkbook()

## =========================================================
## 3) Contrastes vs CN + Tabelas + Volcano
## =========================================================

conds <- c("03dpi","05dpi","10dpi","15dpi","30dpi")

Resumo_tipo_list     <- list()
Resumo_tipo_pct_list <- list()
Resumo_updown_list   <- list()

## guardar DEGs por contraste + união global
DEGS_by_contrast <- list()
ALL_DEGS <- character(0)

for (ct in conds) {
  message("\n=== CONTRASTE: ", ct, " vs CN ===")
  contrast_vec <- c("condition", ct, "CN")
  
  res <- results(dds, contrast = contrast_vec) %>% as.data.frame()
  res$target_id <- rownames(res)
  
  res$target_id_clean <- clean_tx_id(res$target_id)
  
  shr <- do_shrink(dds, contrast_vec)
  if (is(shr, "DESeqResults")) {
    res$log2FC_shrunk <- shr$log2FoldChange
  } else if (is.numeric(shr) && length(shr) == nrow(res)) {
    res$log2FC_shrunk <- as.numeric(shr)
  } else {
    res$log2FC_shrunk <- res$log2FoldChange
  }
  res$log2FC_use <- res$log2FC_shrunk
  
  res$Tipo <- tag_tipo(res$target_id)
  
  res_annot <- res %>%
    left_join(tx2gene, by = c("target_id_clean" = "transcript_id")) %>%
    select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
           log2FC_shrunk, log2FC_use,
           target_id, target_id_clean,
           gene_id, gene_name, gene_biotype,
           Tipo)
  
  res_sig <- res_annot %>%
    filter(!is.na(padj) & padj <= padj_thr & abs(log2FC_use) >= lfc_thr)
  
  ## guardar ids de DEGs
  DEGS_by_contrast[[ct]] <- res_sig$target_id
  ALL_DEGS <- unique(c(ALL_DEGS, res_sig$target_id))
  
  ## salvar TSVs
  write.table(res_annot,
              file = file.path(dir_tabs,
                               sprintf("DESeq2_%s_vs_CN_GTF_all.tsv", ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(res_sig,
              file = file.path(dir_tabs,
                               sprintf("DESeq2_%s_vs_CN_GTF_significant.tsv", ct)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  ## volcano
  if (nrow(res_annot) > 0) {
    plot_volcano(
      df           = res_annot,
      titulo       = sprintf("Volcano - %s vs CN (TEST)", ct),
      outfile_base = file.path(dir_vol, sprintf("Volcano_%s_vs_CN_TEST", ct)),
      padj_thr     = padj_thr,
      lfc_thr      = lfc_thr
    )
  }
  
  ## resumos XM/MSTRG/Other
  tab_tipo <- res_annot %>%
    count(Tipo, name = "n") %>%
    mutate(contrast = paste0(ct, "_vs_CN"))
  
  tab_tipo_pct <- tab_tipo %>%
    group_by(contrast) %>%
    mutate(perc = round(100 * n / sum(n), 2)) %>%
    ungroup()
  
  Resumo_tipo_list[[ct]]     <- tab_tipo
  Resumo_tipo_pct_list[[ct]] <- tab_tipo_pct
  
  ## Up/Down significativos
  if (nrow(res_sig) > 0) {
    updown <- res_sig %>%
      mutate(Direction = if_else(log2FC_use > 0, "Up", "Down")) %>%
      count(Direction, name = "n") %>%
      mutate(contrast = paste0(ct, "_vs_CN"))
  } else {
    updown <- tibble(Direction = c("Up","Down"),
                     n         = 0L,
                     contrast  = paste0(ct, "_vs_CN"))
  }
  Resumo_updown_list[[ct]] <- updown
  
  ## abas Excel
  if (has_xlsx) {
    nm_all <- paste0(ct, "_all")
    nm_sig <- paste0(ct, "_signif")
    openxlsx::addWorksheet(wb, nm_all)
    openxlsx::writeData(wb, nm_all, res_annot)
    openxlsx::addWorksheet(wb, nm_sig)
    openxlsx::writeData(wb, nm_sig, res_sig)
  }
}

## juntar resumos
Resumo_XM_MSTRG_por_contraste <- bind_rows(Resumo_tipo_list)
Resumo_XM_MSTRG_percentual    <- bind_rows(Resumo_tipo_pct_list)
Resumo_Significant_UpDown     <- bind_rows(Resumo_updown_list)

write.table(Resumo_XM_MSTRG_por_contraste,
            file = file.path(dir_tabs, "Resumo_XM_MSTRG_por_contraste.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Resumo_XM_MSTRG_percentual,
            file = file.path(dir_tabs, "Resumo_XM_MSTRG_por_contraste_percentual.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Resumo_Significant_UpDown,
            file = file.path(dir_tabs, "Resumo_Significant_UpDown_por_contraste.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (has_xlsx) {
  openxlsx::addWorksheet(wb, "Resumo_XM_MSTRG")
  openxlsx::writeData(wb, "Resumo_XM_MSTRG", Resumo_XM_MSTRG_por_contraste)
  openxlsx::addWorksheet(wb, "Resumo_XM_MSTRG_%")
  openxlsx::writeData(wb, "Resumo_XM_MSTRG_%", Resumo_XM_MSTRG_percentual)
  openxlsx::addWorksheet(wb, "Resumo_UpDown")
  openxlsx::writeData(wb, "Resumo_UpDown", Resumo_Significant_UpDown)
  
  xlsx_file <- file.path(out_root, "DESeq2_TEST_Complete.xlsx")
  openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  message("Excel salvo em: ", xlsx_file)
}

## =========================================================
## 4) PCA global (DEGs + corte por variância) — SVG + PDF
## =========================================================
message("\n>> PCA global (VST) usando DEGs + corte por variância ...")

if (length(ALL_DEGS) >= 2) {
  degs_ok <- intersect(ALL_DEGS, rownames(vsd))
  if (length(degs_ok) >= 2) {
    mat_deg <- assay(vsd)[degs_ok, , drop = FALSE]
    vv <- matrixStats::rowVars(mat_deg)
    if (length(vv) > PCA_MAX_GENES_VAR) {
      ord <- order(vv, decreasing = TRUE)
      degs_ok <- degs_ok[ord][seq_len(PCA_MAX_GENES_VAR)]
    }
    vsd_pca <- vsd[degs_ok, ]
  } else {
    vsd_pca <- vsd
  }
} else {
  message("(!) Poucos DEGs globais; PCA global será com todos os genes filtrados (rowSums>=2).")
  vsd_pca <- vsd
}

pca_data <- plotPCA(vsd_pca, intgroup = "condition",
                    ntop = min(nrow(vsd_pca), nrow(vsd_pca)),
                    returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% var")) +
  ylab(paste0("PC2: ", percentVar[2], "% var")) +
  theme_minimal(base_size = 14) +
  ggtitle("PCA global - TEST (DEGs + var)")

ggsave(file.path(dir_pca_g, "PCA_Global_DEGs_Var_TEST.svg"),
       p_pca, width = 7, height = 6, device = "svg")
ggsave(file.path(dir_pca_g, "PCA_Global_DEGs_Var_TEST.pdf"),
       p_pca, width = 7, height = 6)

## =========================================================
## 5) PCA por condição vs CN (DEGs do contraste + var)
## =========================================================
message(">> PCA por condição vs CN (DEGs + var) ...")
for (ct in conds) {
  message("   * ", ct)
  try(plot_pca_contraste(ct, data, coldata, dir_pca_c,
                         degs_ct = DEGS_by_contrast[[ct]]),
      silent = TRUE)
}

## =========================================================
## 6) HEATMAPS (com TRIPLICATAS)
##     - GLOBAL: união de DEGs (todos os contrastes) em TODAS as amostras
##     - POR CONTRASTE: apenas CN + condição (com triplicatas)
##     Saída: SVG + PDF
## =========================================================
message(">> Heatmaps — GLOBAL (ALL_DEGS) + por contraste (CN vs ct) ...")

dir_heatmap_global <- mk(dir_hm, "Global_ALL_DEGs")
dir_heatmap_degs   <- mk(dir_hm, "DEGs_por_contraste")

palette <- hcl.colors(256, palette = "Teal")

## 6.1 Heatmap GLOBAL com ALL_DEGS (todas as amostras, triplicatas incluídas)
if (length(ALL_DEGS) >= 2) {
  degs_ok <- intersect(ALL_DEGS, rownames(vsd))
  if (length(degs_ok) >= 2) {
    mat_global <- assay(vsd)[degs_ok, , drop = FALSE]
    
    ## limitar por variância se ficar gigante
    if (!is.null(HEATMAP_MAX_GENES_VAR) && nrow(mat_global) > HEATMAP_MAX_GENES_VAR) {
      vv <- matrixStats::rowVars(as.matrix(mat_global))
      ord <- order(vv, decreasing = TRUE)
      mat_global <- mat_global[ord[seq_len(HEATMAP_MAX_GENES_VAR)], , drop = FALSE]
    }
    
    ## Z-score por gene
    mat_global_z <- t(scale(t(mat_global)))
    
    ann <- data.frame(condition = coldata$condition)
    rownames(ann) <- rownames(coldata)
    
    ## SVG
    svg(file.path(dir_heatmap_global, "Heatmap_Global_ALL_DEGs_TEST.svg"),
        width = 11, height = 8.5)
    pheatmap(
      mat_global_z,
      annotation_col = ann,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      color = palette,
      border_color = NA,
      fontsize_col = 9,
      main = "TEST — Heatmap GLOBAL (ALL_DEGs; triplicatas)"
    )
    dev.off()
    
    ## PDF
    pdf(file.path(dir_heatmap_global, "Heatmap_Global_ALL_DEGs_TEST.pdf"),
        width = 11, height = 8.5)
    pheatmap(
      mat_global_z,
      annotation_col = ann,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      color = palette,
      border_color = NA,
      fontsize_col = 9,
      main = "TEST — Heatmap GLOBAL (ALL_DEGs; triplicatas)"
    )
    dev.off()
    
  } else {
    message("(!) ALL_DEGS não intersecta com vsd (degs_ok < 2).")
  }
} else {
  message("(!) Poucos DEGs globais para heatmap global (ALL_DEGS < 2).")
}

## 6.2 Função heatmap por contraste (CN + ct, triplicatas) — SVG + PDF
make_deg_heatmap <- function(ct, dds, vsd, coldata,
                             alpha = padj_thr,
                             lfc_thresh = 0,
                             max_genes_var = HEATMAP_MAX_GENES_VAR,
                             outdir = dir_heatmap_degs) {
  
  contrast_vec <- c("condition", ct, "CN")
  res <- results(dds, contrast = contrast_vec)
  
  sel <- which(!is.na(res$padj) & res$padj < alpha & abs(res$log2FoldChange) > lfc_thresh)
  degs <- rownames(res)[sel]
  
  if (length(degs) < 2) {
    message("   (!) Poucos DEGs em ", ct, " vs CN (n=", length(degs), "). Pulando heatmap.")
    return(invisible(NULL))
  }
  
  keep <- coldata$condition %in% c("CN", ct)
  keep_samples <- rownames(coldata)[keep]
  
  mat <- assay(vsd)[degs, keep_samples, drop = FALSE]
  
  ## limitar por variância se ficar gigante
  if (!is.null(max_genes_var) && nrow(mat) > max_genes_var) {
    vv <- matrixStats::rowVars(as.matrix(mat))
    ord <- order(vv, decreasing = TRUE)
    mat <- mat[ord[seq_len(max_genes_var)], , drop = FALSE]
  }
  
  mat_z <- t(scale(t(mat)))
  
  ann <- data.frame(condition = coldata$condition)
  rownames(ann) <- rownames(coldata)
  ann <- ann[keep_samples, , drop = FALSE]
  
  out_svg <- file.path(outdir, paste0("Heatmap_DEGs_", ct, "_vs_CN_TEST.svg"))
  out_pdf <- file.path(outdir, paste0("Heatmap_DEGs_", ct, "_vs_CN_TEST.pdf"))
  
  ## SVG
  svg(out_svg, width = 9, height = 9)
  pheatmap(
    mat_z,
    annotation_col = ann,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    color = palette,
    border_color = NA,
    fontsize_col = 10,
    main = paste0("TEST — DEGs ", ct, " vs CN (triplicatas)")
  )
  dev.off()
  
  ## PDF
  pdf(out_pdf, width = 9, height = 9)
  pheatmap(
    mat_z,
    annotation_col = ann,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    color = palette,
    border_color = NA,
    fontsize_col = 10,
    main = paste0("TEST — DEGs ", ct, " vs CN (triplicatas)")
  )
  dev.off()
  
  message(" Heatmap salvo: ", out_svg)
}

message(">> Heatmaps por contraste (CN vs ct) ...")
for (ct in conds) {
  message("   * ", ct, " vs CN")
  try(make_deg_heatmap(ct, dds, vsd, coldata,
                       alpha = padj_thr,
                       lfc_thresh = 0,
                       max_genes_var = HEATMAP_MAX_GENES_VAR,
                       outdir = dir_heatmap_degs),
      silent = TRUE)
}

## =========================================================
## 7) Gráfico-resumo (barras Up/Down por contraste) — SVG + PDF
## =========================================================
try({
  df_bar <- Resumo_Significant_UpDown %>%
    mutate(contrast = factor(contrast, levels = paste0(conds, "_vs_CN")))
  
  p_bar <- ggplot(df_bar, aes(x = contrast, y = n, fill = Direction)) +
    geom_col(position = "stack") +
    theme_minimal(base_size = 12) +
    labs(x = "Contraste",
         y = "Nº DEGs",
         title = "DEGs significativos (padj ≤ 0.05 & |LFC| ≥ 1) — TEST") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(dir_sum, "Barplot_DEGs_UpDown_TEST.svg"),
         p_bar, width = 8.5, height = 5.5, device = "svg")
  ggsave(file.path(dir_sum, "Barplot_DEGs_UpDown_TEST.pdf"),
         p_bar, width = 8.5, height = 5.5)
}, silent = TRUE)

## =========================================================
## 8) Encerrar
## =========================================================
capture.output(sessionInfo(),
               file = file.path(dir_logs, "R_sessionInfo_TEST.txt"))
message("\n Pipeline TEST concluído!")
message("Saídas em: ", out_root)

sink(type = "message"); close(zz)
