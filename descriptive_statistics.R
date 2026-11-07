########Descriptive Statistics########

######Load packages
pacman::p_load(tidyverse, gtsummary)

###### Table 1
dem_clean %>%
  mutate(race = fct_relevel(race, "Unknown", after = Inf)) %>%
  tbl_summary(include = c(biological_sex, race, age),
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 2,
              label = list(biological_sex = "Sex", race = "Race", age = "Age")) %>%
              modify_header(label ~ "**Demographic**")
            
  