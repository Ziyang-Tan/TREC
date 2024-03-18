library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(EnhancedVolcano)

data <- read_csv("data/processed/TREC_KREC_baby_combine_all.csv")
subtype_names <- colnames(data)[24:35]
protein_names <- colnames(data)[39:126]

data %>%
    filter(!is.na(Neutrophils)) %>%
    dim() # check the number of samples we have both data

ggplot(data, aes(x = week, y = `LOG10 sjTREC`)) +
    geom_point()

lapply(c("birth_group", "group", "sex", "mode_delivery", "birth_hospital"), function(x) {
    ggplot(data %>% filter(!is.na(CD8T)), aes_string(x = "CD8T", y = "`LOG10 sjTREC`", color = x)) +
        geom_point() +
        theme_bw() +
        theme(legend.position = "bottom", )
}) %>%
    ggarrange(plotlist = .) %>%
    ggexport(filename = "figures/correlation analysis/log10TREC QC.pdf")

lapply(c("birth_group", "group", "sex", "mode_delivery", "birth_hospital"), function(x) {
    ggplot(data %>% filter(!is.na(CD8T)), aes_string(x = "CD8T", y = "`LOG10 KREC`", color = x)) +
        geom_point() +
        theme_bw() +
        theme(legend.position = "bottom", )
}) %>%
    ggarrange(plotlist = .) %>%
    ggexport(filename = "figures/correlation analysis/log10KREC QC.pdf")

# TREC analysis
# CyTOF subpopulation
res_subtype <- lapply(subtype_names, function(name) {
    cor.test(as.formula(paste0("~ `LOG10 sjTREC` + ", name)), data, use = "complete.obs", method = "spearman")
})
df_subtype <- data.frame(
    r = sapply(res_subtype, function(x) x$estimate),
    p = sapply(res_subtype, function(x) x$p.value),
    subtype = subtype_names
)
EnhancedVolcano(df_subtype,
    x = "r", y = "p",
    lab = df_subtype$subtype, xlim = c(-0.75, 0.75),
    pCutoff = 1, FCcutoff = 0,
    drawConnectors = TRUE,
    max.overlaps = 20,
    title = "TREC and subpopulation frequencies (CyTOF)",
    subtitle = ""
) + xlab("Spearman correlation")
ggsave("figures/correlation analysis/TREC correlation CyTOF.pdf")

lapply(subtype_names, function(sub_name) {
    ggplot(data %>% filter(!is.na(sub_name)), aes_string(x = sub_name, y = "`LOG10 sjTREC`")) +
        geom_point() +
        theme_bw()
}) %>%
    ggarrange(plotlist = .) %>%
    ggexport(filename = "figures/correlation analysis/log10TREC subtypes.pdf")

# Olink protein profile
na_count <- data %>%
    select(all_of(protein_names)) %>%
    is.na() %>%
    colSums()
hist(na_count / dim(data)[1])
protein_names_qc <- protein_names[na_count / dim(data)[1] <= 0.9]
res_protein <- lapply(protein_names_qc, function(name) {
    cor.test(as.formula(paste0("~ `LOG10 sjTREC` + `", name, "`")), data, use = "complete.obs", method = "spearman")
})
df_protein <- data.frame(
    r = sapply(res_protein, function(x) x$estimate),
    p = sapply(res_protein, function(x) x$p.value),
    name = protein_names_qc
)
EnhancedVolcano(df_protein,
    x = "r", y = "p",
    lab = df_protein$name, xlim = c(-0.75, 0.75),
    pCutoff = 1, FCcutoff = 0,
    drawConnectors = TRUE,
    max.overlaps = 20,
    title = "TREC and protein profiles (Olink)",
    subtitle = "",
    arrowheads = FALSE
) + xlab("Spearman correlation")
ggsave("figures/correlation analysis/TREC correlation Olink.pdf", width = 12, height = 12)

select_proteins <- c("TRANCE", "TNFB", "CXCL5", "IL7", "IL8", "IL6", "`MMP-10`", "IL10")
lapply(select_proteins, function(sub_name) {
    ggplot(data %>% filter(!is.na(sub_name)), aes_string(x = sub_name, y = "`LOG10 sjTREC`")) +
        geom_point() +
        theme_bw()
}) %>%
    ggarrange(plotlist = .) %>%
    ggexport(filename = "figures/correlation analysis/log10TREC proteins.pdf")

ggplot(data %>% filter(!is.na(`TRANCE`)), aes(x = `TRANCE`, y = `LOG10 sjTREC`)) +
    geom_point() +
    theme_bw()
ggsave("figures/correlation analysis/log10TREC TRANCE.pdf")

write_csv(df_protein, "data/processed/TREC_correlation_proteins.csv")
write_csv(df_subtype, "data/processed/TREC_correlation_celltypes.csv")

# KREC analysis
res_subtype <- lapply(subtype_names, function(name) {
    cor.test(as.formula(paste0("~ `LOG10 KREC` + ", name)), data, use = "complete.obs", method = "spearman")
})
df_subtype <- data.frame(
    r = sapply(res_subtype, function(x) x$estimate),
    p = sapply(res_subtype, function(x) x$p.value),
    subtype = subtype_names
)
EnhancedVolcano(df_subtype,
    x = "r", y = "p",
    lab = df_subtype$subtype, xlim = c(-0.75, 0.75),
    pCutoff = 1, FCcutoff = 0,
    drawConnectors = TRUE,
    max.overlaps = 20,
    title = "KREC and subpopulation frequencies (CyTOF)",
    subtitle = ""
) + xlab("Spearman correlation")
ggsave("figures/correlation analysis/KREC correlation CyTOF.pdf")

res_protein <- lapply(protein_names_qc, function(name) {
    cor.test(as.formula(paste0("~ `LOG10 KREC` + `", name, "`")), data, use = "complete.obs", method = "spearman")
})
df_protein <- data.frame(
    r = sapply(res_protein, function(x) x$estimate),
    p = sapply(res_protein, function(x) x$p.value),
    name = protein_names_qc
)
EnhancedVolcano(df_protein,
    x = "r", y = "p",
    lab = df_protein$name, xlim = c(-0.4, 0.4),
    pCutoff = 1, FCcutoff = 0,
    drawConnectors = TRUE,
    max.overlaps = 20,
    title = "KREC and protein profiles (Olink)",
    subtitle = "",
    arrowheads = FALSE
) + xlab("Spearman correlation")
ggsave("figures/correlation analysis/KREC correlation Olink.pdf", width = 12, height = 12)

write_csv(df_protein, "data/processed/KREC_correlation_proteins.csv")
write_csv(df_subtype, "data/processed/KREC_correlation_celltypes.csv")
