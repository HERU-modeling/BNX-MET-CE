rm(list = ls()) # to clean the workspace

library(dplyr) # to manipulate data
library(reshape2) # to transform data
library(ggplot2) # for nice looking plots
# library(ggbreak)
library(tidyverse)
# library(data.table)
library(tidyr)
library(RColorBrewer)

## Load DSA output files
# QALYs
load(file = "outputs/dsa/owsa/outcomes_owsa_itt_ps.RData")
load(file = "outputs/dsa/owsa/outcomes_owsa_pp.RData")

# Load data sets
df_owsa_itt_ps <- l_owsa_itt_ps$df_owsa_itt_ps
df_owsa_pp <- l_owsa_pp$df_owsa_pp

# DSA labels
df_dsa_ly_labels_itt_ps <- read.csv(file = "data/dsa/owsa/tornado_labels_ly_itt_ps.csv", header = TRUE)
df_dsa_ly_labels_pp <- read.csv(file = "data/dsa/owsa/tornado_labels_ly_pp.csv", header = TRUE)

# Subset by mean
# Deterministic
n_base_itt_ps <- l_owsa_itt_ps$n_base_itt_ps
n_base_pp <- l_owsa_pp$n_base_pp

## Combine data frames
# QALYs
# ITT-PS
df_owsa_itt_ps <- as_tibble(df_owsa_itt_ps) %>%
  mutate(base = n_base_itt_ps) %>%
  mutate(diff = ifelse(abs(Upper - Lower) > 0, abs(Upper - Lower), abs(base - Upper)))

# PP
df_owsa_pp <- as_tibble(df_owsa_pp) %>%
  mutate(base = n_base_pp) %>%
  mutate(diff = ifelse(abs(Upper - Lower) > 0, abs(Upper - Lower), abs(base - Upper)))

#########################
#### Tornado diagram ####
#########################
### QALYs ###
# ITT-PS
v_order_parameters <- df_owsa_itt_ps %>%
  arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>%
  unlist() %>%
  levels()

# width of columns in plot (value between 0 and 1)
width <- 0.6
# get data frame in shape for ggplot and geom_rect
df.2 <- df_owsa_itt_ps %>%
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key = "type", value = "output.value", Lower:Upper) %>%
  # just reordering columns
  select(var_name, type, output.value, diff, base) %>%
  # create the columns for geom_rect
  mutate(
    var_name = factor(var_name, levels = v_order_parameters),
    ymin = pmin(output.value, base),
    ymax = pmax(output.value, base),
    xmin = as.numeric(var_name) - width / 2,
    xmax = as.numeric(var_name) + width / 2
  )
# Add value labels for SA ranges
data_merge2 <- inner_join(df.2, df_dsa_ly_labels_itt_ps, by = c("var_name", "type")) # Applying inner_join() function

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
p_tornado_itt_ps <- ggplot() +
  geom_rect(
    data = data_merge2,
    aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)
  ) +
  geom_text(data = data_merge2, aes(y = output.value, x = (xmax + xmin) / 2, label = dsa_itt_ps_u), hjust = 1) +
  geom_text(data = data_merge2, aes(y = output.value, x = (xmax + xmin) / 2, label = dsa_itt_ps_l), hjust = 0) +
  theme_bw() +
  scale_fill_manual(values = c(
    "Upper" = "midnightblue",
    "Lower" = "slategray2"
  )) +
  theme(
    axis.title.y = element_blank(), legend.position = "none",
    legend.title = element_blank()
  ) +
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(
    breaks = c(seq_along(v_order_parameters)),
    labels = v_order_parameters
  ) +
  xlab("Parameter") +
  ylab("Incremental life years (BNX vs. methadone)") +
  ylim(-3500, 0) +
  coord_flip()

# PP
v_order_parameters <- df_owsa_pp %>%
  arrange(diff) %>%
  mutate(var_name = factor(x = var_name, levels = var_name)) %>%
  select(var_name) %>%
  unlist() %>%
  levels()

# width of columns in plot (value between 0 and 1)
width <- 0.6
# get data frame in shape for ggplot and geom_rect
df.2 <- df_owsa_pp %>%
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key = "type", value = "output.value", Lower:Upper) %>%
  # just reordering columns
  select(var_name, type, output.value, diff, base) %>%
  # create the columns for geom_rect
  mutate(
    var_name = factor(var_name, levels = v_order_parameters),
    ymin = pmin(output.value, base),
    ymax = pmax(output.value, base),
    xmin = as.numeric(var_name) - width / 2,
    xmax = as.numeric(var_name) + width / 2
  )
# Add value labels for SA ranges
data_merge2 <- inner_join(df.2, df_dsa_ly_labels_pp, by = c("var_name", "type")) # Applying inner_join() function

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)
p_tornado_pp <- ggplot() +
  geom_rect(
    data = data_merge2,
    aes(ymax = ymax, ymin = ymin, xmax = xmax, xmin = xmin, fill = type)
  ) +
  geom_text(data = data_merge2, aes(y = output.value, x = (xmax + xmin) / 2, label = dsa_pp_u), hjust = 1) +
  geom_text(data = data_merge2, aes(y = output.value, x = (xmax + xmin) / 2, label = dsa_pp_l), hjust = 0) +
  theme_bw() +
  scale_fill_manual(values = c(
    "Upper" = "midnightblue",
    "Lower" = "slategray2"
  )) +
  theme(
    axis.title.y = element_blank(), legend.position = "none",
    legend.title = element_blank()
  ) +
  geom_hline(yintercept = df.2$base) +
  scale_x_continuous(
    breaks = c(seq_along(v_order_parameters)),
    labels = v_order_parameters
  ) +
  xlab("Parameter") +
  ylab("Incremental life years (BNX vs. methadone)") +
  ylim(-3500, 0) +
  coord_flip()

# Output plots
ggsave("Plots/dsa/owsa/tornado_itt_ps.png", p_tornado_itt_ps, height = 5, width = 7, dpi = 320)
ggsave("Plots/dsa/owsa/tornado_pp.png", p_tornado_pp, height = 5, width = 7, dpi = 320)
