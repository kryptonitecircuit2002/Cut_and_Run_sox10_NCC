library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

load_file <- function(file, variant, tp) {
  df <- read_tsv(
    file,
    col_names = FALSE,
    na = c(".")
  )
  colnames(df) <- c("chr","start","end","b_chr","b_start","b_end","b_name","b_score", "b_strand", "signalValue", "pvalue",
                    "qvalue", "peak")
  df <- df %>%
    mutate(
      peak_id = paste(chr, start, end, sep="_"),
      Variant = variant,
      TP = tp,
      active = ifelse(is.na(b_chr), 0, 1)
    ) %>%
    dplyr::select(peak_id, TP, Variant, active, signalValue)
  
  return(df)
}

files <- list(
  list(".../h2_peaks_18hpf.bed","H2A.Z","16hpf"),
  list(".../h2_peaks_24hpf.bed","H2A.Z","24hpf"),
  list(".../h2_peaks_48hpf.bed","H2A.Z","48hpf"),
  list(".../h3_peaks_18hpf.bed","H3.3","16hpf"),
  list(".../h3_peaks_24hpf.bed","H3.3","24hpf"),
  list(".../h3_peaks_48hpf.bed","H3.3","48hpf")
)

variant_list <- lapply(files, function(x) load_file(x[[1]], x[[2]], x[[3]]))
variant_all <- bind_rows(variant_list)
variant_all <- variant_all %>%
  group_by(peak_id, TP, Variant) %>%
  summarise(active = max(active), .groups = "drop")
variant_wide <- variant_all %>%
  pivot_wider(
    names_from = Variant,
    values_from = active,
    values_fill = 0
  )

classify_variant <- function(H2A.Z, H3.3) {
  if (H2A.Z == 1 & H3.3 == 0) return("H2A.Z")
  if (H2A.Z == 0 & H3.3 == 1) return("H3.3")
  if (H2A.Z == 1 & H3.3 == 1) return("H2A.Z+H3.3")
  "Unassigned"
}
variant_wide$State <- mapply(classify_variant,
                             variant_wide$H2A.Z,
                             variant_wide$H3.3)
variant_levels <- c("H2A.Z", "H3.3", "H2A.Z+H3.3", "No Variant")

variant_matrix <- variant_wide %>%
  dplyr::select(peak_id, TP, State) %>%
  pivot_wider(
    names_from = TP,
    values_from = State,
    values_fill = "No Variant"
  )

df_alluvial <- variant_matrix %>%
  count(`16hpf`, `24hpf`, `48hpf`, name = "Freq") %>%
  complete(
    `16hpf` = variant_levels,
    `24hpf` = variant_levels,
    `48hpf` = variant_levels,
    fill = list(Freq = 0)
  )


state_colors <- c(
  "H2A.Z" = "#4C72B0",  
  "H3.3" = "#DD8452",  
  "H2A.Z+H3.3" = "#55A868", 
  "Non variant" = "thistle3"
)


ggplot(
  df_alluvial,
  aes(
    axis1 = `16hpf`,
    axis2 = `24hpf`,
    axis3 = `48hpf`,
    y = Freq
  )
) +
  geom_alluvium(aes(fill = `16hpf`), width = 1/6, curve_type = "sigmoid") +
  
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6) +
  
  scale_fill_manual(values = state_colors, drop = FALSE) +
  
  scale_x_discrete(
    limits = c("16hpf", "24hpf", "48hpf"),
    labels = c("16 hpf", "24 hpf", "48 hpf")
  ) +
  
  theme_void() +
  theme(
    text = element_text(family = "Times New Roman", size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    y = "Peaks",
    x = "",
    fill = "Variant State"
  )

