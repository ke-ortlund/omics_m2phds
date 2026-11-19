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
mon_dem_clean %>%
  mutate(race = fct_relevel(race, "Unknown", after = Inf)) %>%
  tbl_summary(by = study_accession, include = c(biological_sex, race, age, value_clean),
              type = list(age = "continuous2", value_clean = "continuous2"),
              statistic = list(
                all_continuous() ~ c("{mean} ({sd})", "{median} ({p25}, {p75})"),
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = list(
                all_continuous() ~ 1,
                all_categorical() ~ c(0, 1)),
              label = list(biological_sex = "Sex", race = "Race", age = "Age")) %>%
              modify_header(label ~ "**Demographic**")

###### S        
flow_mon_clean %>%
  ggplot(aes(x = value_clean)) + 
  geom_histogram(fill = "#1874CD") +
  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
  coord_cartesian(xlim = c(0, 35)) +
  labs(x = "Classical monocytes (%/parent population) ", y = "Frequency")

