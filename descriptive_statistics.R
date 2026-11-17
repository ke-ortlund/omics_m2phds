###############################################################################
# Descriptive Statistics
###############################################################################


###### Load packages
pacman::p_load(tidyverse, gtsummary)

###### Set global theme
theme_set(theme_bw() +
            theme(title = element_text(size = 15, color = "grey30"),
                  axis.title.x = element_text(size = 15, color = "grey30"),
                  axis.title.y = element_text(size = 15, color = "grey30"),
                  axis.text.x = element_text(size = 15, color = "grey30"),
                  axis.text.y = element_text(size = 15, color = "grey30"),
                  panel.grid = element_blank(),
                  strip.text = element_text(size = 15, color = "grey30"),
                  strip.background = element_blank(), 
                  legend.position = "bottom",
                  legend.key.size = unit(3, "cm"),
                  legend.text = element_text(size = 15, color = "grey30"))) 


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

###### S        
flow_mon_clean %>%
  ggplot(aes(x = value_clean)) + 
  geom_histogram() +
  labs(x = "Classical monocytes (%/parent population) ", y = "Frequency")


