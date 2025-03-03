library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

setwd("~/Documents/GitHub/dupcoal/test_MLMSC/")

dup_results <- read.csv('results/fungi_results.csv')
disc_results <- read.csv('results/fungi_results_disco.csv')

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

dupplot <- ggplot(duplications_plot_data, aes(x = lb, y = mean_value, 
                                              color = program, linetype = as.factor(c))) +
  geom_line(lwd=1) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) +
  facet_wrap(~ ld_label, ncol = 2, scales = "free_y") +  # Use the new ld_label column as strings
  labs(x = expression(lambda), y = "Duplications", color = "Program", linetype = "c") +  # Lambda for x-axis label
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

copyplot <- ggplot(gene_copies_plot_data, aes(x = lb, y = mean_value, 
                                  color = program, linetype = as.factor(c))) +
  geom_line(lwd=1) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) +
  facet_wrap(~ ld_label, ncol = 2, scales = "free_y") +  # Use the new ld_label column as strings
  labs(x = expression(lambda), y = "Gene Copies", color = "Program", linetype = "c") +  
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

png("DuplicationsandGeneCopies.png", height=800, width=1000)
grid.arrange(dupplot, copyplot, ncol = 1)
dev.off()

# discordance
head(disc_results)
# summarize data
summary_disc_results <- disc_results %>%
  pivot_longer(cols = c(total_quartets, normalized_quartet_score), 
               names_to = "metric", 
               values_to = "value") %>%
  group_by(program, metric, observed_duplications) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop")

# total quartets
quartets_plotdata <- summary_disc_results %>%
  filter(metric == "total_quartets")

plot_total_quartets <- ggplot(quartets_plotdata, aes(x = observed_duplications, y = mean_value, color = program)) +
  geom_line(lwd=1) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) +
  labs(x = "Duplications", y = "Total Quartets", color = "Program") +  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 18),  # Larger x-axis label
    axis.title.y = element_text(size = 18),  # Larger y-axis label
    legend.text = element_text(size = 16),    # Larger legend text
    legend.title = element_text(size=18),
    axis.text = element_text(size=16)
  )


# total quartets
nquartets_plotdata <- summary_disc_results %>%
  filter(metric == "normalized_quartet_score")

plot_normalized_quartets <- ggplot(nquartets_plotdata, aes(x = observed_duplications, y = mean_value, color = program)) +
  geom_line(lwd=1) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1) +
  labs(x = "Duplications", y = "Normalized Quartet Score", color = "Program") +  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 18),  # Larger x-axis label
    axis.title.y = element_text(size = 18),  # Larger y-axis label
    legend.text = element_text(size = 16),    # Larger legend text
    legend.title = element_text(size=18),
    axis.text = element_text(size=16)
  )
# Output both plots in a single figure
pdf("Discordance_summary.pdf", height=10, width=8)
grid.arrange(plot_total_quartets, plot_normalized_quartets, ncol = 1)
dev.off()
