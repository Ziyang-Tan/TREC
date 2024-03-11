library(dplyr)
library(readr)

raw <- read_csv("data/raw_TRECdata_combined.csv.gz")
cytof_data <- read_csv("data/all_exp_lineageFreq_wideformat_wMeta.csv")
olink_data <- NULL

# fetch meta data from bigquery
raw <- read_delim("data/BilanTELLUSmerged_edit.csv", delim = ";") %>%
    mutate(across(contains("QTY"), .fns = function(x) {
        as.numeric(gsub(" ", "", x))
    }))
baby_id <- raw %>%
    filter(type == "Baby") %>%
    select(study_id) %>%
    distinct() %>%
    unlist(use.names = F)
sql <- paste0("
SELECT mode_delivery, baby.* FROM `cradle-259115.cradle.subjects_delivery` AS Delivery
LEFT JOIN `cradle-259115.cradle.subjects_baby` AS baby
on Delivery.baby_id = baby.subject
WHERE baby.study_id IN (", paste(baby_id, collapse = ","), ")")
tb <- bq_project_query("cradle-259115", sql)
tmp <- bq_table_download(tb)
raw <- raw %>% left_join(tmp, by = "study_id")

# QC for TREC
data <- raw %>%
    filter(
        log10(`sjTREC QTY`) > 1, # QC
        type == "Baby"
    ) %>%
    mutate(
        sex = case_match(sex, 0 ~ "male", 1 ~ "female"),
        child_id = as.character(child_id),
        log10TREC = log10(`sjTREC QTY`),
        mode_delivery = case_match(
            mode_delivery,
            c("Acute C-section", "Planned C-section") ~ "C-section",
            c("Vaginal-instrumental", "Vaginal") ~ "Vaginal"
        )
    )
