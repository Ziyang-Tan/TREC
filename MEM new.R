library(dplyr)
library(readr)
library(ggplot2)
library(lme4)
library(bigrquery)
library(ggpubr)
library(rstatix)

data <- read_csv("data/processed/TREC_baby_combine_all.csv")
df_protein <- read_csv("data/processed/correlation_proteins.csv")
df_subtype <- read_csv("data/processed/correlation_celltypes.csv")

subtype_names <- colnames(data)[22:33]
protein_names <- colnames(data)[37:124]
highcor_protein_names <- df_protein %>%
    arrange(desc(abs(r))) %>%
    select(name) %>%
    slice(1:10)

data_qc <- data %>% filter(log10TREC > 1)

data_logit <- data_qc %>%
    mutate(normlog10TREC = (log10TREC - min(log10TREC)) / (max(log10TREC) - min(log10TREC)))
df <- data_logit %>% select(all_of(highcor_protein_names$name), log10TREC, week, `GENOTYPE SNP rs2204985`, sex, CD4T, CD8T, pDC, Neutrophils)
m <- glm(log10TREC ~ week + `GENOTYPE SNP rs2204985` + sex + CD4T + CD8T + pDC + Neutrophils + TRANCE, data = df, na.action = "na.omit")
m <- glm(log10TREC ~ week + `GENOTYPE SNP rs2204985` + sex + CD4T + CD8T, data = df, na.action = "na.omit")

summary(m)
m$model$fitted <- predict(m, type = "response")
ggplot(m$model) +
    geom_point(aes(week, log10TREC)) +
    geom_point(aes(week, fitted, color = "red")) +
    theme_bw() +
    theme(panel.background = element_blank())
performance::r2(m)
