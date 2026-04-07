## =========================================================
##  Immune analysis — TESTÍCULO 
##  Autor: Walter + ChatGPT
##  Saídas: SVG + PDF
## =========================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
  
  ## enriquecimento
  library(clusterProfiler)
  library(org.Mm.eg.db)  # (mesma estratégia do seu TEST anterior)
  
  ## rótulos bonitos no volcano
  library(ggrepel)
})

## ----------------------------
## (A) CAMINHOS TESTÍCULO
## ----------------------------
out_root    <- "/Users/walterfranco/Library/CloudStorage/OneDrive-Pessoal(2)/Doutorado iec/Experimentos/Analise_seq/analisados/test/test_22_12"
dir_tabs    <- file.path(out_root, "tabelas")
matriz_file <- "/Users/walterfranco/Library/CloudStorage/OneDrive-Pessoal(2)/Doutorado iec/Experimentos/Analise_seq/analisados/test/matriz_kallisto.isoform.counts.matrix"

## onde salvar tudo da parte imune (TESTÍCULO)
mk <- function(...) { d <- file.path(...); dir.create(d, showWarnings=FALSE, recursive=TRUE); d }
dir_immune <- mk(out_root, "immune_analysis")

## ----------------------------
## (B) CONDIÇÕES TESTÍCULO
## ----------------------------
conds <- c("03dpi","05dpi","10dpi","15dpi","30dpi")

## ----------------------------
## Funções
## ----------------------------
get_condition <- function(s) {
  s <- toupper(s)
  if (grepl("^CN", s)) return("CN")
  m <- regmatches(s, regexpr("^\\d{2}DPI", s))
  if (length(m) == 1 && nchar(m) > 0) return(gsub("DPI","dpi",tolower(m)))
  "CN"
}

read_counts_matrix <- function(path){
  if (!file.exists(path)) stop("Arquivo não encontrado: ", path)
  d <- read.table(path, header=TRUE, row.names=1, sep="\t",
                  check.names=FALSE, quote="", comment.char="")
  as.matrix(d)
}

## =========================================================
## 1) Recriar dds/vsd/coldata do TESTÍCULO
## =========================================================
data <- read_counts_matrix(matriz_file)
data <- round(data)
data <- data[rowSums(data) >= 2, , drop=FALSE]

samples   <- colnames(data)
condition <- vapply(samples, get_condition, character(1))

## garante fatores consistentes
condition <- factor(condition, levels = c(conds, "CN"))
coldata   <- data.frame(condition = condition, row.names = samples)
stopifnot(ncol(data) == nrow(coldata))

dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

## =========================================================
## 2) Enriquecimento GO (BP) usando DEGs (por contraste)
##    -> constrói immune_gene_symbols
## =========================================================

## 2.1) juntar genes DEGs (gene_name) de TODOS os contrastes
all_sig_paths <- file.path(dir_tabs, paste0("DESeq2_", conds, "_vs_CN_GTF_significant.tsv"))
names(all_sig_paths) <- conds

sig_all <- lapply(all_sig_paths, function(p) {
  if (!file.exists(p)) return(NULL)
  readr::read_tsv(p, show_col_types = FALSE)
})
sig_all <- sig_all[!vapply(sig_all, is.null, logical(1))]

if (length(sig_all) == 0) stop("Nenhuma tabela significant encontrada em: ", dir_tabs)

all_deg_symbols <- bind_rows(sig_all) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::pull(gene_name) %>%
  unique()

message("Total gene_name (DEGs anotados) para enriquecimento: ", length(all_deg_symbols))

## 2.2) SYMBOL -> ENTREZ (mouse, como no seu workflow)
gene_df <- bitr(all_deg_symbols,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Mm.eg.db)

if (is.null(gene_df) || nrow(gene_df) == 0) {
  stop("Falha no bitr(): nenhum mapeamento SYMBOL->ENTREZID retornado. Verifique símbolos e OrgDb.")
}

## 2.3) GO BP
ego_bp <- enrichGO(gene          = unique(gene_df$ENTREZID),
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

if (is.null(ego_bp) || nrow(ego_bp@result) == 0) stop("enrichGO não retornou termos.")

## 2.4) filtrar termos “imunes” por palavras-chave
immune_terms <- ego_bp@result %>%
  dplyr::filter(grepl("immune|immun|cytokine|interferon|inflamm|antigen|defense|leukocyte|chemokine|antiviral",
                      Description, ignore.case = TRUE))

message("Termos GO (BP) com palavras-chave imunes: ", nrow(immune_terms))

immune_gene_symbols <- unique(unlist(strsplit(immune_terms$geneID, "/")))
immune_gene_symbols <- immune_gene_symbols[immune_gene_symbols != ""]
message("Genes (símbolos) presentes nesses termos imunes: ", length(immune_gene_symbols))

## salva lista para referência
writeLines(immune_gene_symbols, file.path(dir_immune, "immune_gene_symbols_TESTICULO.txt"))

## =========================================================
## 3) Heatmap imune GLOBAL (todas as condições, triplicatas)
## =========================================================
dir_global <- mk(dir_immune, "GLOBAL")
dir_global_hm <- mk(dir_global, "heatmap")

## pegar TODOS os IDs (target_id) que correspondem a genes imunes em QUALQUER contraste
immune_target_ids_global <- bind_rows(sig_all) %>%
  dplyr::filter(!is.na(gene_name), gene_name %in% immune_gene_symbols) %>%
  dplyr::select(target_id, gene_name) %>%
  dplyr::distinct()

immune_ids_global <- base::intersect(immune_target_ids_global$target_id, rownames(vsd))

if (length(immune_ids_global) >= 2) {
  
  mat <- assay(vsd)[immune_ids_global, , drop=FALSE]
  
  ## mapear gene_name para linhas
  gene_map <- immune_target_ids_global
  new_rn <- gene_map$gene_name[match(rownames(mat), gene_map$target_id)]
  rownames(mat) <- new_rn
  
  ## limpar NAs e duplicados
  mat <- mat[!is.na(rownames(mat)), , drop=FALSE]
  mat <- mat[!duplicated(rownames(mat)), , drop=FALSE]
  
  ## z-score por gene
  mat_z <- t(scale(t(mat)))
  
  ann <- data.frame(condition = coldata$condition)
  rownames(ann) <- rownames(coldata)
  
  pal <- colorRampPalette(rev(brewer.pal(11, "Blues")))(256)
  
  ## SVG (vetorial) — base R
   svg(file.path(dir_global_hm, "Heatmap_immune_GLOBAL_TESTICULO11.svg"),
      width = 11, height = 11)
  pheatmap(mat_z,
           annotation_col = ann,
           show_rownames  = TRUE,
           #fontsize_row   = ifelse(nrow(mat_z) > 50, 5, 7),
           fontsize_row   = 15,
           fontsize_col   = 10,
           color          = pal,
           border_color   = NA,
           main           = "TESTÍCULO - Immune genes (GLOBAL, triplicatas)")
  dev.off()
  
  ## PDF (backup)
  pdf(file.path(dir_global_hm, "Heatmap_immune_GLOBAL_TESTICULO.pdf"),
      width = 11, height = 11)
  pheatmap(mat_z,
           annotation_col = ann,
           show_rownames  = TRUE,
           fontsize_row   = ifelse(nrow(mat_z) > 50, 5, 7),
           fontsize_col   = 10,
           color          = pal,
           border_color   = NA,
           main           = "TESTÍCULO - Immune genes (GLOBAL, triplicatas)")
  dev.off()
  
} else {
  message("(!) Poucos immune IDs para heatmap GLOBAL (n=", length(immune_ids_global), ")")
}

## =========================================================
## 4) Por contraste: tabela + volcano rotulado + heatmap (CN+ct)
## =========================================================
for (ct in conds) {
  
  message("\n=== IMUNE: ", ct, " vs CN ===")
  
  dir_ct     <- mk(dir_immune, paste0(ct, "_vs_CN"))
  dir_ct_tab <- mk(dir_ct, "tabelas")
  dir_ct_vol <- mk(dir_ct, "volcano")
  dir_ct_hm  <- mk(dir_ct, "heatmap")
  
  sig_path <- file.path(dir_tabs, paste0("DESeq2_", ct, "_vs_CN_GTF_significant.tsv"))
  if (!file.exists(sig_path)) {
    message("(!) Não achei: ", sig_path, " — pulando ", ct)
    next
  }
  
  sig_ct <- readr::read_tsv(sig_path, show_col_types = FALSE)
  
  immune_ct <- sig_ct %>%
    dplyr::filter(!is.na(gene_name)) %>%
    dplyr::filter(gene_name %in% immune_gene_symbols)
  
  ## salvar tabela
  readr::write_tsv(immune_ct, file.path(dir_ct_tab, paste0("DEGs_immune_", ct, "_vs_CN_TESTICULO.tsv")))
  
  ## volcano com rótulos (top)
  if (nrow(immune_ct) >= 2) {
    
    immune_ct2 <- immune_ct %>%
      dplyr::mutate(Signif = ifelse(!is.na(padj) & padj <= 0.05, "padj<=0.05", "NS"))
    
    label_df <- immune_ct2 %>%
      dplyr::filter(
        padj <= 0.01 |
          abs(log2FC_use) >= quantile(abs(log2FC_use), 0.90, na.rm = TRUE)
      ) %>%
      dplyr::distinct(gene_name, .keep_all = TRUE)
    
    p <- ggplot(immune_ct2, aes(x = log2FC_use, y = -log10(pvalue), color = Signif)) +
      geom_point(alpha = 0.7, size = 1.6) +
      ggrepel::geom_text_repel(data = label_df, aes(label = gene_name),
                               size = 3, max.overlaps = 30) +
      theme_minimal(base_size = 12) +
      labs(title = paste0("Volcano (immune DEGs) - ", ct, " vs CN (TESTÍCULO)"),
           x = "log2FC (shrink)", y = "-log10(p-value)")
    
    ## SVG (vetorial) — sem svglite
    ggsave(file.path(dir_ct_vol, paste0("Volcano_immune_", ct, "_vs_CN_TESTICULO.svg")),
           p, width = 7, height = 6, device = "svg")
    
    ## PDF (backup)
    ggsave(file.path(dir_ct_vol, paste0("Volcano_immune_", ct, "_vs_CN_TESTICULO.pdf")),
           p, width = 7, height = 6)
    
  } else {
    message("(!) Poucos immune DEGs para volcano em ", ct, " (n=", nrow(immune_ct), ")")
  }
  
  ## heatmap (CN + ct, triplicatas) com gene_name nas linhas
  immune_ids <- unique(immune_ct$target_id)
  immune_ids <- base::intersect(immune_ids, rownames(vsd))
  
  if (length(immune_ids) >= 2) {
    
    keep_samples <- rownames(coldata)[coldata$condition %in% c("CN", ct)]
    mat <- assay(vsd)[immune_ids, keep_samples, drop=FALSE]
    
    gene_map <- immune_ct %>%
      dplyr::select(target_id, gene_name) %>%
      dplyr::distinct()
    
    rownames(mat) <- gene_map$gene_name[match(rownames(mat), gene_map$target_id)]
    mat <- mat[!is.na(rownames(mat)), , drop=FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop=FALSE]
    
    ## z-score por gene
    mat_z <- t(scale(t(mat)))
    
    ann <- data.frame(condition = coldata$condition)
    rownames(ann) <- rownames(coldata)
    ann <- ann[keep_samples, , drop=FALSE]
    
    pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)
    
    ## SVG (vetorial)
    svg(file.path(dir_ct_hm, paste0("Heatmap_immune_", ct, "_vs_CN_TESTICULO.svg")),
        width = 10, height = 10)
    pheatmap(mat_z,
             annotation_col = ann,
             show_rownames  = TRUE,
             fontsize_row   = ifelse(nrow(mat_z) > 50, 5, 7),
             fontsize_col   = 10,
             color          = pal,
             border_color   = NA,
             main           = paste0("Immune DEGs - ", ct, " vs CN (TESTÍCULO, triplicatas)"))
    dev.off()
    
    ## PDF (backup)
    pdf(file.path(dir_ct_hm, paste0("Heatmap_immune_", ct, "_vs_CN_TESTICULO.pdf")),
        width = 10, height = 10)
    pheatmap(mat_z,
             annotation_col = ann,
             show_rownames  = TRUE,
             fontsize_row   = ifelse(nrow(mat_z) > 50, 5, 7),
             fontsize_col   = 10,
             color          = pal,
             border_color   = NA,
             main           = paste0("Immune DEGs - ", ct, " vs CN (TESTÍCULO, triplicatas)"))
    dev.off()
    
  } else {
    message("(!) Poucos immune IDs p/ heatmap em ", ct, " (n=", length(immune_ids), ")")
  }
}

message("\n✅ Immune analysis TESTÍCULO finalizada! Saída em: ", dir_immune)
