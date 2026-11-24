###############################################################################
# Descriptive Statistics
###############################################################################


###### Load packages
pacman::p_load(tidyverse, gtsummary)

###### Table 1
table_1 <- mon_dem_clean %>%
  mutate(race = fct_recode(factor(race), Missing = "Unknown"),
         vacc_any = if_else(vacc_any %in% 1, "Participants with at least 1 vaccination", "Missing"),
         across(c(race, vacc_any), ~ fct_relevel(.x, "Missing", after = Inf))) %>%
  tbl_summary(by = study_accession, include = c(biological_sex, race, age, vacc_any, value_clean),
              type = list(age = "continuous2", value_clean = "continuous2"),
              statistic = list(
                all_continuous() ~ c("{mean} ({sd})", "{median} ({p25}, {p75})"),
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = list(
                all_continuous() ~ 1,
                all_categorical() ~ c(0, 1)),
              label = list(biological_sex = "Sex", 
                           race = "Race", 
                           age = "Age", 
                           vacc_any = "Vaccine history", 
                           value_clean = "Classical Monocyte"),
              missing_text = "Missing") %>%
              modify_header(label ~ "") %>%
              bold_labels() %>%
              add_overall() %>%
              modify_footnote(all_stat_cols() ~ NA)


# Save
table_1 %>%
  as_gt() %>%
  gt::gtsave(filename = "table_1.docx", path = output_folder)
  
###### Figure 1 
figure_1 <- mon_dem_clean %>%
  ggplot(aes(x = study_accession, y = value_clean, fill = study_accession)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.6) +
  scale_y_continuous(breaks = seq(0, 40, by = 5), limits = c(0, NA)) +
  scale_x_discrete(labels = c("2024_UGA" = "UGA 2024", "SDY2867" = "SDY 2867")) +
  labs(x = "Study",
       y = "Classical monocytes (percent of parent population, Day 7)") +
  theme(legend.position = "none")

# Save
ggsave(filename = file.path(output_folder, "figure_1.png"),
       plot = figure_1, width = 8, height = 7, dpi = 300, bg = "white")

###### Normality check

# Histogram 
histogram <- mon_dem_clean %>%
  ggplot(aes(x = value_clean)) +
  geom_histogram(aes(),
                 color = "white", binwidth = 2, alpha = 0.8) +
  geom_density(aes(y = after_stat(count) * 2),
               linewidth = 1) +
  scale_fill_manual(values = study_cols) +
  scale_color_manual(values = line_cols) +
  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
  coord_cartesian(xlim = c(0, 35)) +
 # facet_wrap(~ study_accession,
 #            labeller = as_labeller(c("2024_UGA" = "UGA 2024",
#                                    "SDY2867" = "SDY2867"))) +
  labs(x = "Classical monocytes (percent of parent population, Day 7)",
       y = "Count") +
  theme(legend.position = "none")


ggsave(filename = file.path(output_folder, "histogram.png"),
       plot = histogram, width = 8, height = 7, dpi = 300, bg = "white")

# Shapiro-Wilk Test
shapiro.test(mon_dem_clean$value_clean)

# Q-Q Plot
qqnorm(mon_dem_clean$value_clean,
       main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles",
       pch = 19,
       col = "blue")

qqline(mon_dem_clean$value_clean,
       col = "red",
       lwd = 2)

