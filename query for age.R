library(dplyr)
library(readr)
library(bigrquery)


dat <- read_delim("data/BilanTELLUSmerged_edit.csv", delim = ";")

mother_id <- dat %>%
  filter(type == "Mother") %>%
  select(study_id) %>%
  unlist(use.names = F)
sql <- paste0("SELECT age, study_id FROM `cradle-259115.cradle.subjects_mother` WHERE study_id IN (", paste(mother_id, collapse = ","), ")")
tb <- bq_project_query("cradle-259115", sql)
mother_data <- bq_table_download(tb)

father_id <- dat %>%
  filter(type == "Father") %>%
  select(study_id) %>%
  unlist(use.names = F)
sql <- paste0("SELECT age, study_id FROM `cradle-259115.cradle.subjects_father` WHERE study_id IN (", paste(father_id, collapse = ","), ")")
tb <- bq_project_query("cradle-259115", sql)
father_data <- bq_table_download(tb)

tmp <- dat %>%
  left_join(rbind(mother_data, father_data), by = "study_id") %>%
  mutate(age = case_when(
    type %in% c("Mother", "Father") & age == 0 ~ NA,
    TRUE ~ age
  )) %>%
  rename(adult_age = age)

write_csv(tmp, "data/BilanTELLUSmerged_parent_age.csv")


# all samples
mother_sample_id <- dat %>%
  filter(type == "Mother") %>%
  mutate(ref_id = paste0('"', referral_id, '"')) %>% # add "" for the SQL query below.
  select(ref_id) %>%
  unlist(use.names = F)
sql <- paste0("SELECT * FROM `cradle-259115.cradle.all_samples_AronEdit` WHERE sample_id IN (", paste(mother_sample_id, collapse = ","), ")")
tb <- bq_project_query("cradle-259115", sql)
df <- bq_table_download(tb)

df %>% filter(sample_id == "2600048978")
