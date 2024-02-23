library(dplyr)
library(readr)
library(ggplot2)
library(lme4)
library(bigrquery)
library(ggpubr)
library(rstatix)

# raw <- read_delim('data/BilanTELLUSmerged_edit.csv', delim = ';') %>%
#   mutate(across(contains('QTY'), .fns= function(x) {as.numeric(gsub(' ', '', x))}))

# # fetch info from bigquery
# # sex, age, mode of delivery for babies. 
# baby_id <- raw %>% filter(type == 'Baby') %>% select(study_id) %>% distinct() %>% unlist(use.names = F)
# sql <- paste0("
# SELECT mode_delivery, baby.* FROM `cradle-259115.cradle.subjects_delivery` AS Delivery 
# LEFT JOIN `cradle-259115.cradle.subjects_baby` AS baby
# on Delivery.baby_id = baby.subject
# WHERE baby.study_id IN (", paste(baby_id, collapse = ','), ')')
# 
# tb <- bq_project_query('cradle-259115', sql)
# tmp <- bq_table_download(tb)
# raw <- raw %>% left_join(tmp, by='study_id')
# 
# write_csv(raw, 'data/raw_TRECdata_combined.csv.gz')

raw <- read_csv('data/raw_TRECdata_combined.csv.gz')
cytof_data <- read_csv('data/all_exp_lineageFreq_wideformat_wMeta.csv')

ggplot(raw, aes(x=week, y=log10(`sjTREC QTY`))) + geom_point() + theme_bw()

data <- raw %>% 
  filter(log10(`sjTREC QTY`) >1, # QC
         type == 'Baby') %>%
  mutate(sex = case_match(sex, 0 ~ 'male', 1 ~ 'female'),
         child_id = as.character(child_id),
         log10TREC = log10(`sjTREC QTY`),
         mode_delivery = case_match(mode_delivery, 
                                    c('Acute C-section', 'Planned C-section') ~ 'C-section',
                                    c('Vaginal-instrumental', 'Vaginal') ~ 'Vaginal'
         ))
adata <- data %>%
  mutate(sample_id = as.character(referral_id)) %>%
  left_join(cytof_data %>% 
              select(all_of(c('sample_id', 'B.cells', 'Basophils', 'CD4T', 'CD8T', 'gdT', 'Monocytes', 
                              'Neutrophils', 'NK.cells', 'pDC', 'mDC', 'Tregs', 'Eosinophils', 
                              'Platelets', 'MSC'
              ))), by='sample_id') %>% 
  filter(type == 'Baby') %>%
  mutate(birth_group = case_when(
    week < 5 ~ 'birth',
    week >=5 & week < 29 ~ 'early',
    week >= 29 ~ 'late'
  ))

ggplot(adata %>% filter(!is.na(CD4T)), aes(x=CD4T, y=log10TREC)) + geom_point() + theme_bw()

# try to reproduce some figures
ggplot(data, aes(x=week, y=log10TREC)) + geom_point() + theme_bw()
# ggplot(data, aes(x=week, y=log10TREC, color=`GENOTYPE SNP rs2204985`)) + geom_point() + geom_smooth()
# ggplot(adata, aes(x=birth_group, y=log10TREC)) + geom_boxplot() + geom_jitter()
# ggplot(adata, aes(x=birth_group, y=log10TREC, color = `GENOTYPE SNP rs2204985`)) + geom_boxplot()

# MEM
m <- lmer(log10TREC ~ week + `GENOTYPE SNP rs2204985` + sex + mode_delivery + (1|child_id), data)
tdat <- data.frame(predicted=predict(m, type = 'response'), residual = residuals(m), child_id = data$child_id)
ggplot(tdat,aes(x=predicted,y=residual)) + geom_point() + geom_hline(yintercept=0, lty=3)
qqnorm(tdat$residual)
ggplot(cbind(data, tdat$predicted), aes(x=week, y=log10TREC))+ geom_point() + geom_point(aes(x=week, y=`tdat$predicted`, color='red'))
rstatix::Anova(m)
performance::r2(m)

# m <- lm(log10TREC ~ week + `GENOTYPE SNP rs2204985` + sex + mode_delivery, data)

ggplot(data, aes(x=week, y=log10TREC)) + geom_point() + theme_bw()

adataLogit <- adata %>% 
  mutate(normlog10TREC = (log10TREC - min(log10TREC)) / (max(log10TREC) - min(log10TREC)))

m <- glm(normlog10TREC ~ week + `GENOTYPE SNP rs2204985` + sex + mode_delivery, family = binomial(link = "probit"), data = adataLogit)
summary(m)
m$model$fitted <- predict(m, type = "response")
ggplot(m$model) + 
  geom_point(aes(week, normlog10TREC)) +
  geom_point(aes(week, fitted, color='red')) +
  theme_bw() + 
  theme(panel.background = element_blank())
performance::r2(m)

m <- gam::gam(log10TREC ~ week + `GENOTYPE SNP rs2204985` + sex + mode_delivery, data = adata)
m <- glm(log10TREC ~ week + `GENOTYPE SNP rs2204985` + sex + mode_delivery, data=adata)

d.stats <- adata %>% t_test(log10TREC ~ sex, p.adjust.method = 'fdr') %>% add_xy_position(x = 'sex')

ggplot(adata, aes(x=sex, y=log10TREC)) + geom_boxplot() + geom_jitter() +stat_pvalue_manual(d.stats, tip.length = 0.01, label='p') + theme_bw()
ggsave('figures/TREC by sex.pdf', width = 5, height = 5)

g1 <- ggplot(adata, aes(x=week, y=log10TREC, color=sex)) + geom_point() + geom_smooth()
d.stats <- adata %>% group_by(birth_group) %>% t_test(log10TREC ~ sex, p.adjust.method = 'fdr') %>% add_xy_position(x = 'birth_group')
g2 <- ggplot(adata, aes(x=birth_group, y=log10TREC, color=sex)) + geom_boxplot()+stat_pvalue_manual(d.stats, tip.length = 0.01, label='p') + theme_bw()
ggarrange(g1, g2)
ggsave('figures/TREC by sex and age.pdf', width = 10, height = 5)

ggplot(adata %>% filter(!is.na(CD8T)), aes(x=CD8T, y=log10TREC)) + geom_point() + theme_bw()

# KREC
data <- raw %>%
  filter(#log10(`sjTREC QTY`) >1, # QC
    !is.na(`sjKREC QTY`),
    type == 'Baby') %>%
  mutate(sex = case_match(sex, 0 ~ 'male', 1 ~ 'female'),
         child_id = as.character(child_id),
         log10KREC = log10(`sjKREC QTY` + 1),
         mode_delivery = case_match(mode_delivery,
                                    c('Acute C-section', 'Planned C-section') ~ 'C-section',
                                    c('Vaginal-instrumental', 'Vaginal') ~ 'Vaginal'
         ))
adata <- data %>%
  mutate(sample_id = as.character(referral_id)) %>%
  left_join(cytof_data %>% 
              select(all_of(c('sample_id', 'B.cells', 'Basophils', 'CD4T', 'CD8T', 'gdT', 'Monocytes', 
                              'Neutrophils', 'NK.cells', 'pDC', 'mDC', 'Tregs', 'Eosinophils', 
                              'Platelets', 'MSC'
              ))), by='sample_id') %>% 
  filter(type == 'Baby') %>%
  mutate(birth_group = case_when(
    week < 5 ~ 'birth',
    week >=5 & week < 29 ~ 'early',
    week >= 29 ~ 'late'
  ))
ggplot(adata, aes(x=week, y=log10KREC)) + geom_point()

ggplot(adata, aes(x=week, y=log10KREC)) + geom_point() + geom_smooth()

m <- glm(log10KREC ~ week + `GENOTYPE SNP rs2204985` + sex + mode_delivery, data = adata)
tdat <- data.frame(predicted=predict(m, type = 'response'), residual = residuals(m), child_id = data$child_id)
ggplot(tdat,aes(x=predicted,y=residual)) + geom_point() + geom_hline(yintercept=0, lty=3)
qqnorm(tdat$residual)
ggplot(cbind(data, tdat$predicted), aes(x=week, y=log10KREC))+ geom_point() + geom_point(aes(x=week, y=`tdat$predicted`, color='red'))
rstatix::Anova(m)
performance::r2(m)

m <- lmer(log10KREC ~ week + `GENOTYPE SNP rs2204985` + sex + mode_delivery + (1|child_id), data = adata)
tdat <- data.frame(predicted=predict(m, type = 'response'), residual = residuals(m), child_id = data$child_id)
ggplot(tdat,aes(x=predicted,y=residual)) + geom_point() + geom_hline(yintercept=0, lty=3)
qqnorm(tdat$residual)
ggplot(cbind(data, tdat$predicted), aes(x=week, y=log10KREC))+ geom_point() + geom_point(aes(x=week, y=`tdat$predicted`, color='red'))
rstatix::Anova(m)
performance::r2(m)

ggplot(adata %>% filter(!is.na(B.cells)), aes(x=B.cells, y=log10KREC)) + geom_point() + theme_bw()




