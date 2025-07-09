### PREPROCESS AND PLOT FLUORESCENCE MICROSCOPY DATA ###

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

conditions <- c("F", "Fcvx", "Fm", "Fmcvx")
t <- seq(0, 49.75, length.out = 200)

all_data <- list()

# Extract and preprocess data
for (cond in conditions) {
  data <- read_csv(paste0(cond, ".csv"))
  nROI <- (ncol(data) - 1) / 3
  
  roi_dfs <- list()
  for (i in 1:nROI) {
    mean_col <- paste0("Mean", i)
    area_col <- paste0("Area", i)
    ch_col <- paste0("Ch", i)
    
    # Normalize signal
    norm_signal <- 1000 * data[[mean_col]] / data[[area_col]]
    
    # Split by channel
    ch <- data[[ch_col]]
    for (ch_type in c(1, 3)) {
      ch_signal <- norm_signal[ch == ch_type]
      if (length(ch_signal) == length(t)) {
        roi_dfs[[length(roi_dfs) + 1]] <- data.frame(
          time = t,
          signal = ch_signal,
          ROI = paste0("ROI", i),
          Channel = ifelse(ch_type == 1, "Aux/IAA", "ARF"),
          Condition = cond
        )
      }
    }
  }
  
  all_data[[cond]] <- bind_rows(roi_dfs)
}

# Combine all data
plot_data <- bind_rows(all_data)

# Summarize to mean and SD per time point
summary_data <- plot_data %>%
  group_by(Condition, Channel, time) %>%
  summarise(
    mean_signal = mean(signal, na.rm = TRUE),
    sd_signal = sd(signal, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot
ggplot(summary_data, aes(x = time, y = mean_signal, color = Channel, fill = Channel)) +
  geom_ribbon(aes(ymin = mean_signal - sd_signal, ymax = mean_signal + sd_signal), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  facet_wrap(~Condition, scales = "free_y") +
  scale_color_manual(values = c("ARF" = "#00BFC4", "Aux/IAA" = "#F8766D")) +
  scale_fill_manual(values = c("ARF" = "#00BFC4", "Aux/IAA" = "#F8766D")) +
  labs(title = "Signal over Time by Condition and Channel",
       x = "Time (min)", y = "Normalized Signal") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))



# Filter data for F and Fcvx
comparison_data <- summary_data %>%
  filter(Condition %in% c("F", "Fcvx"))

# Plot with ggplot
ggplot(comparison_data, aes(x = time, y = mean_signal,
                            color = Channel, linetype = Condition, fill = Channel)) +
  geom_ribbon(aes(ymin = mean_signal - sd_signal, ymax = mean_signal + sd_signal), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Aux/IAA" = "#F8766D", "ARF" = "#00BFC4")) +
  scale_fill_manual(values = c("Aux/IAA" = "#F8766D", "ARF" = "#00BFC4")) +
  scale_linetype_manual(values = c("F" = "solid", "Fcvx" = "dashed"), labels = c("F" = "- cvxIAA", "Fcvx" = "+ cvxIAA")) +
  labs(title = "ACwt", x = "Time (h)", y = "Normalized Signal") +
  theme_minimal()

# Filter data for F and Fcvx
comparison_data <- summary_data %>%
  filter(Condition %in% c("Fm", "Fmcvx"))

# Plot with ggplot
ggplot(comparison_data, aes(x = time, y = mean_signal,
                            color = Channel, linetype = Condition, fill = Channel)) +
  geom_ribbon(aes(ymin = mean_signal - sd_signal, ymax = mean_signal + sd_signal), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Aux/IAA" = "#F8766D", "ARF" = "#00BFC4")) +
  scale_fill_manual(values = c("Aux/IAA" = "#F8766D", "ARF" = "#00BFC4")) +
  scale_linetype_manual(values = c("Fm" = "solid", "Fmcvx" = "dashed"), labels = c("Fm" = "- cvxIAA", "Fmcvx" = "+ cvxIAA")) +
  labs(title = "ACm3", x = "Time (h)", y = "Normalized Signal") +
  theme_minimal()

peak_data <- plot_data %>%
  group_by(Condition, Channel, ROI) %>%
  summarise(PeakValue = max(signal, na.rm = TRUE), .groups = "drop")

# F vs Fcvx
F_vs_Fcvx <- peak_data %>%
  filter(Condition %in% c("F", "Fcvx"))

# By channel
wilcox_F_aux <- wilcox.test(PeakValue ~ Condition, data = filter(F_vs_Fcvx, Channel == "Aux/IAA"))
wilcox_F_arf <- wilcox.test(PeakValue ~ Condition, data = filter(F_vs_Fcvx, Channel == "ARF"))

# Fm vs Fmcvx
Fm_vs_Fmcvx <- peak_data %>%
  filter(Condition %in% c("Fm", "Fmcvx"))

# By channel
wilcox_Fm_aux <- wilcox.test(PeakValue ~ Condition, data = filter(Fm_vs_Fmcvx, Channel == "Aux/IAA"))
wilcox_Fm_arf <- wilcox.test(PeakValue ~ Condition, data = filter(Fm_vs_Fmcvx, Channel == "ARF"))

# Print results
wilcox_F_aux$p.value
wilcox_F_arf$p.value
wilcox_Fm_aux$p.value
wilcox_Fm_arf$p.value

# Plot all conditions together
ggplot(peak_data, aes(x = Condition, y = PeakValue, fill = Condition)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  facet_wrap(~Channel, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#B3CDE3", "Fcvx" = "#6497B1", "Fm" = "#FBB4AE", "Fmcvx" = "#F768A1")) +
  scale_x_discrete(labels = c("F" = "ACwt \n- cvxIAA", 
                            "Fcvx" = "ACwt \n+ cvxIAA", 
                            "Fm" = "ACm3 \n- cvxIAA", 
                            "Fmcvx" = "ACm3 \n+ cvxIAA")) +
  labs(title = "Peak Fluorescence Across All Samples",
       x = "Sample",
       y = "Peak Normalized Signal") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

#############
#Extra plots#
#############

# Pivot wider to get Aux/IAA and ARF signals side by side for each ROI, time, condition
ratio_data <- plot_data %>%
  pivot_wider(names_from = Channel, values_from = signal) %>%
  filter(!is.na(`Aux/IAA`) & !is.na(ARF)) %>%
  mutate(Ratio = `Aux/IAA` / ARF)

summary_ratio <- ratio_data %>%
  group_by(Condition, time) %>%
  summarise(
    mean_ratio = mean(Ratio, na.rm = TRUE),
    sd_ratio = sd(Ratio, na.rm = TRUE),
    .groups = 'drop'
  )

ggplot(summary_ratio, aes(x = time, y = mean_ratio, color = Condition, fill = Condition)) +
  geom_ribbon(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Aux/IAA to ARF Signal Ratio Over Time by Condition",
       x = "Time (min)", y = "Mean Aux/IAA / ARF Ratio") +
  theme_minimal()

peak_ratio <- ratio_data %>%
  group_by(Condition, ROI) %>%
  summarise(PeakRatio = max(Ratio, na.rm = TRUE), .groups = 'drop')

ggplot(peak_ratio, aes(x = Condition, y = PeakRatio, fill = Condition)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Peak Aux/IAA to ARF Signal Ratio by Condition",
       x = "Condition", y = "Peak Aux/IAA / ARF Ratio") +
  theme_minimal() +
  theme(legend.position = "none")

