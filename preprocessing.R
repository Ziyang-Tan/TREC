library(dplyr)
library(readr)
library(bigrquery)

cytof_data <- read_csv("data/all_exp_lineageFreq_wideformat_wMeta.csv")
olink_data <- read_delim("data/Olink_TRECsubsampling.csv", delim = ";")

# fetch meta data from bigquery
# raw <- read_delim("data/BilanTELLUSmerged_edit.csv", delim = ";") %>%
#     mutate(across(contains("QTY"), .fns = function(x) {
#         as.numeric(gsub(" ", "", x))
#     }))
raw <- read_csv2("data/TELLUS_all_data_corrected_032024.csv") %>%
    mutate(
        `LOG10 sjTREC` = as.numeric(sub("\\,", "\\.", `LOG10 sjTREC`)),
        `LOG10 CjInt` = as.numeric(sub("\\,", "\\.", `LOG10 CjInt`)),
        `LOG10 KREC` = as.numeric(sub("\\,", "\\.", `LOG10 KREC`))
    )
baby_id <- raw %>%
    filter(type == "Baby") %>%
    select(study_id) %>%
    distinct() %>%
    unlist(use.names = FALSE)
sql <- paste0("
SELECT mode_delivery, baby.* FROM `cradle-259115.cradle.subjects_delivery` AS Delivery
LEFT JOIN `cradle-259115.cradle.subjects_baby` AS baby
on Delivery.baby_id = baby.subject
WHERE baby.study_id IN (", paste(baby_id, collapse = ","), ")")
tb <- bq_project_query("cradle-259115", sql)
tmp <- bq_table_download(tb)
raw <- raw %>% left_join(tmp, by = "study_id")

data <- raw %>%
    filter(
        type == "Baby"
    ) %>%
    mutate(
        sex = case_match(sex, 0 ~ "male", 1 ~ "female"),
        child_id = as.character(child_id),
        mode_delivery = case_match(
            mode_delivery,
            c("Acute C-section", "Planned C-section") ~ "C-section",
            c("Vaginal-instrumental", "Vaginal") ~ "Vaginal"
        )
    )

olink <- olink_data %>%
    filter(MissingFreq <= 0.5) %>%
    tidyr::pivot_wider(
        names_from = "Assay", values_from = "NPX",
        id_cols = "SampleID", values_fn = ~ mean(.x, na.rm = TRUE)
    ) %>%
    rename(sample_id = SampleID) %>%
    mutate(sample_id = as.character(sample_id))

adata <- data %>%
    mutate(sample_id = as.character(referral_id)) %>%
    left_join(cytof_data %>%
        select(all_of(c(
            "sample_id", "B.cells", "Basophils", "CD4T", "CD8T", "gdT", "Monocytes",
            "Neutrophils", "NK.cells", "pDC", "mDC", "Tregs", "Eosinophils",
            "Platelets", "MSC", "baby_age"
        ))), by = "sample_id") %>%
    left_join(olink, by = "sample_id") %>%
    filter(type == "Baby") %>%
    mutate(birth_group = case_when(
        week < 5 ~ "birth",
        week >= 5 & week < 29 ~ "early",
        week >= 29 ~ "late"
    ))

adata %>%
    filter(!is.na(Neutrophils)) %>%
    dim() # check the number of samples we have both data
write_csv(adata, "data/processed/TREC_KREC_baby_combine_all.csv")
