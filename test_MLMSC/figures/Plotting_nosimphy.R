library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(plotrix)

setwd("~/Documents/GitHub/dupcoal/test_MLMSC/")

dup_results <- read.csv('results/fungi_results.csv')
dup_results <- dup_results[dup_results$program!="simphy",]
# summarize data
summary_dup_results <- dup_results %>%
  pivot_longer(cols = c(observed_duplications, gene_copies), 
               names_to = "metric", 
               values_to = "value") %>%
  group_by(lb, program, c, ld, metric) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop")

# duplications
duplications_plot_data <- summary_dup_results %>%
  filter(metric == "observed_duplications") %>%
  mutate(
    c = ifelse(c == 9, "bl / 9", as.character(c)),
    ld_label = recode(ld, `0.5` = "μ = 0.5 * λ", `lb` = "μ = λ")  # Use strings for facet labels
  )

duplications_plot_data <- duplications_plot_data %>%
  mutate(panel_label = ifelse(ld_label == unique(ld_label)[1], "A", "B"))

dupplot <- ggplot(duplications_plot_data, aes(x = lb, y = mean_value, 
                                              color = program, linetype = as.factor(c))) +
  geom_line(lwd=1) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) +
  facet_wrap(~ ld_label, ncol = 2, scales = "free_y") +  # Use the new ld_label column as strings
  labs(x = expression(lambda), y = "Duplications", color = "Program", linetype = "c") +  # Lambda for x-axis label
  geom_text(aes(x = -Inf, y = Inf, label = panel_label),
            hjust = -0.3, vjust = 1.3, size = 8,
            inherit.aes = FALSE) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 18),  # Adjust facet label size
    axis.title.x = element_text(size = 18),  # Larger x-axis label
    axis.title.y = element_text(size = 18),  # Larger y-axis label
    legend.text = element_text(size = 16),    # Larger legend text
    legend.title = element_text(size=18),
    axis.text = element_text(size=16)
  )

dupplot
# gene copies
gene_copies_plot_data <- summary_dup_results %>%
  filter(metric == "gene_copies") %>%
  mutate(
    c = ifelse(c == 9, "bl / 9", as.character(c)),
    ld_label = recode(ld, `0.5` = "μ = 0.5 * λ", `lb` = "μ = λ")  # Use strings for facet labels
  )

gene_copies_plot_data <- gene_copies_plot_data %>%
  mutate(panel_label = ifelse(ld_label == unique(ld_label)[1], "C", "D"))

copyplot <- ggplot(gene_copies_plot_data, aes(x = lb, y = mean_value, 
                                  color = program, linetype = as.factor(c))) +
  geom_line(lwd=1) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) +
  facet_wrap(~ ld_label, ncol = 2, scales = "free_y") +  # Use the new ld_label column as strings
  labs(x = expression(lambda), y = "Gene Copies", color = "Program", linetype = "c") +  
  geom_text(aes(x = -Inf, y = Inf, label = panel_label),
            hjust = -0.3, vjust = 1.3, size = 8,
            inherit.aes = FALSE) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 18),  # Adjust facet label size
    axis.title.x = element_text(size = 18),  # Larger x-axis label
    axis.title.y = element_text(size = 18),  # Larger y-axis label
    legend.text = element_text(size = 16),    # Larger legend text
    legend.title = element_text(size=18),
    axis.text = element_text(size=16)
  )
copyplot

png("./figures/DuplicationsandGeneCopies_nosimphy.png", height=800, width=1000)
grid.arrange(dupplot, copyplot, ncol = 1)
dev.off()


pdf("./figures/DuplicationsandGeneCopies_nosimphy.pdf", height=8, width=10)
grid.arrange(dupplot, copyplot, ncol = 1)
dev.off()

