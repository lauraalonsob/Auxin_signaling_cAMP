### PREPROCESS, DETREND, PLOT AND COMPARE MICROFLUIDICS DATA ###

# Load packages
library(dplyr)
library(ggplot2)
library(pracma)
library(readr)
library(reshape2)

# Function to extract ROI data by channel
extract_channel_data <- function(df, channel = 3) {
  roi_count <- (ncol(df) - 1) / 3
  values <- list()
  for (i in 1:roi_count) {
    area_col <- 1 + (i - 1) * 3 + 1
    mean_col <- area_col + 1
    ch_col <- mean_col + 1
    roi_data <- df[, c(area_col, mean_col, ch_col)]
    colnames(roi_data) <- c("Area", "Mean", "Ch")
    roi_data <- roi_data %>% filter(Ch == channel)
    values[[paste0("ROI", i)]] <- roi_data$Mean / roi_data$Area
  }
  as.data.frame(values)
}

# Variable settings
file_path <- "microF.csv"        # filename
t <- 10                          # Time between frames (e.g. 10 seconds or minutes)
cutoff_start <- 150     # Frames to remove from the beginning    (150 for microF.csv)
cutoff_end   <- 0     # Frames to remove from the end      (230 for microFm.csv)
channel_to_use <- 1              # 1 for mCherry, 3 for GFP

### LOAD & PREPROCESS DATA ###

# Read raw data
raw <- read_csv(file_path, col_names = FALSE, skip = 1)

# Extract data for the selected channel
df <- extract_channel_data(raw, channel = channel_to_use)

# Apply cutoffs
n_total <- nrow(df)
if (cutoff_start > 0) {
  df <- df[-(1:cutoff_start), ]
}
if (cutoff_end > 0) {
  df <- df[1:(nrow(df) - cutoff_end), ]
}

time <- seq(1, nrow(df) * t, t)
df$time <- time

### DETREND & NORMALIZE EACH ROI ###


for (roi in colnames(df)[!colnames(df) %in% "time"]) {
  signal <- df[[roi]]
  signal <- detrend(signal)
  signal <- savgol(signal, 15)
  signal <- signal + abs(min(signal))
  df[[roi]] <- signal / max(signal)
}

# Reshape for plotting
df_m <- melt(df, id.vars = "time")

# Plot each ROI
ggplot(df_m, aes(x = time, y = value, color = variable)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y") +
  theme_classic() +
  labs(title = "ACwt   Aux/IAA signal",
       subtitle = "Detrended and Normalized ROI Traces",
       x = "Time (min)", y = "Normalized Intensity") +
  theme(legend.position = "none")


###  COMPUTE FFT PERIOD SPECTRA ###

# Individual FFT period spectrum for each ROI
fft_all <- list()
n_points <- nrow(df)
roi_names <- colnames(df)[!colnames(df) %in% "time"]

for (roi in roi_names) {
  signal <- df[[roi]]
  fft_res <- Mod(fft(signal))
  n <- length(fft_res)
  
  period_df <- tibble(
    ROI = roi,
    freq_index = 1:(n / 2),
    period = (n * t) / (1:(n / 2)),
    magnitude = fft_res[2:(n / 2 + 1)]
  )
  fft_all[[roi]] <- period_df
}

# Combine all into one dataframe
fft_combined <- bind_rows(fft_all)

# Plot individual FFT period spectra
ggplot(fft_combined, aes(x = period, y = magnitude)) +
  geom_line(color = "black", size = 0.8) +
  facet_wrap(~ ROI, scales = "free_y") +
  scale_x_continuous(trans = "log10") +
  labs(title = "ACwt   Aux/IAA signal",
       subtitle = "FFT Period Spectrum by ROI",
       x = "Period (min)",
       y = "Magnitude") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Compute mean magnitude across ROIs per period
fft_mean <- fft_combined %>%
  group_by(freq_index, period) %>%
  summarise(mean_magnitude = mean(magnitude, na.rm = TRUE), .groups = "drop")

# Plot mean FFT spectrum
ggplot(fft_mean, aes(x = period, y = mean_magnitude)) +
  geom_line(color = "blue", size = 1) +
  scale_x_continuous(trans = "log10") +
  labs(title = "ACwt   Aux/IAA signal",
       subtitle = "Mean FFT Period Spectrum",
       x = "Period (min)",
       y = "Mean Magnitude") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Store mean FFT spectrum
fft_meanFI=fft_mean

### COMPARISON BETWEEN CONDITIONS ###

# Add condition names
fft_meanFA$Condition <- "ACwt  ARF"
fft_meanFI$Condition <- "ACwt  Aux/IAA"
fft_meanFmA$Condition <- "ACm3  ARF"
fft_meanFmI$Condition <- "ACm3  Aux/IAA"

# Combine all into one dataframe
fft_all_means <- bind_rows(fft_meanFA, fft_meanFI, fft_meanFmA, fft_meanFmI)

# Plot
ggplot(fft_all_means, aes(x = period, y = mean_magnitude, color = Condition)) +
  geom_vline(xintercept = 240, linetype = "dashed", color = "black", size = 0.8) +
  geom_line(size = 1) +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = c(
    "ACwt  ARF" = "cyan2",
    "ACwt  Aux/IAA" = "red2",
    "ACm3  ARF" = "cyan3",
    "ACm3  Aux/IAA" = "red3"
  )) +
  labs(
    title = "Mean FFT Period Spectra by Condition",
    x = "Period (min)",
    y = "Mean Magnitude"
  ) +
  facet_wrap(~ Condition, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "none",  # optional if you don't want a legend
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )