# Data analysis for Lab Report 2.


# 0. Load libraries & import data ----------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)

setwd("./Report2_RT-qPCR/")

# RNA isolations
yields = read.csv(
  './data/RNA_yield_purity_for_plot.csv',
  stringsAsFactors = T
  )

# Primer standardization
primers = read.csv('./data/primer_std_values.csv')

# Cycle threshold
ct_values = read.csv(
  './data/RT_qPCR_final_CT_values.csv',
  stringsAsFactors = T
  )


# 1. Boxplot of tissue masses --------------------------------------------------

mass_plot = yields %>%
  ggplot() +
  geom_boxplot(aes(x = Genotype, y = grams)) +
  facet_grid(~Day) +
  scale_y_continuous(n.breaks = 10) +
  labs(y = 'Tissue mass (g)' )

# Export
ggsave(
  mass_plot,
  './figures/fig1.png',
  units = 'in',
  width = 4.5,
  height = 2.5
  )


# 2. Boxplots: RNA concentrations and yields -----------------------------------

concentration_plot = yields %>%
  mutate(uguL = concentration*0.001) %>%
  ggplot() +
  geom_boxplot(aes(x = Genotype, y = uguL),) +
  facet_grid(~Day) +
  scale_y_continuous(n.breaks = 10) +
  labs(
    y = 'Concentration (μg/μL)',
    title = ''
    )

yield_plot = yields %>%
  # ng/uL * 50 uL elution * (1ng=0.001ug) / initial mass g
  mutate(ug = concentration*50*0.001) %>%
  mutate(yield = ug/(grams/0.1)) %>%
  ggplot() +
  geom_boxplot(aes(x = Genotype, y = yield),) +
  geom_hline(yintercept = 25, linetype = 2) +
  facet_grid(~Day) +
  scale_y_continuous(n.breaks = 8) +
  labs(
    y = 'Yield (μg per 0.1 g tissue)',
    title = ''
    )

concentration_yield_comb = plot_grid(
  concentration_plot,
  yield_plot,
  labels = c('A)', 'B)')
  )

# Export
ggsave(
  concentration_yield_comb,
  './figures/fig2.png',
  units = 'in',
  width = 7.5,
  height = 3.5
)


# 3. Boxplot: Nanodrop absorbance values ---------------------------------------

protein_purity_plot = yields %>%
  ggplot() +
  geom_boxplot(
    aes(x = Genotype, y = protein_purity)
    ) +
  facet_grid(~Day) +
  scale_y_continuous(n.breaks = 10) +
  geom_hline(
    yintercept = 1.8,
    linetype = 2
    ) +
  labs(
    y = 'Absorbance (260/280 nm)',
    title = ''
    )

carb_purity_plot = yields %>%
  ggplot() +
  geom_boxplot(
    aes(x = Genotype, y = carb_purity)
  ) +
  facet_grid(~Day) +
  scale_y_continuous(n.breaks = 7) +
  geom_hline(
    yintercept = 1.8,
    linetype = 2
    ) +
  labs(
    y = 'Absorbance (260/230 nm)',
    title = ''
    )

protein_carb_comb = plot_grid(
  protein_purity_plot,
  carb_purity_plot,
  labels = c('A)', 'B)')
)

# Export
ggsave(
  protein_carb_comb,
  './figures/fig3.png',
  units = 'in',
  width = 7.5,
  height = 3.5
)


# 4. Primer standardization curve ----------------------------------------------

primer_std_curve = primers %>%
  ggplot(aes(x = log10_ng, y = CT)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(
    x = expression('log'[10]*' ng of cDNA'),
    y = 'Cycle Threshold'
    ) +
  stat_regline_equation(
    label.x = -0.5,
    label.y = 33,
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))
    )

# Export
ggsave(
  primer_std_curve,
  './figures/fig5.png',
  units = 'in',
  width = 7,
  height = 3
)


# 5. Sample sizes of valid RT-qPCR CT values, after removing outliers ----------

sample_size_plot = ct_values %>%
  ggplot() +
  geom_bar(
    aes(x = genotype, fill = gene),
    position = 'dodge',
    alpha = 0.9
    ) +
  facet_grid(~age) +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_manual(
    values = c('#5B7553', '#C4AF9A', '#81523F')
    ) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    x = 'Genotype',
    y = 'Sample counts'
    )

# Export
ggsave(
  sample_size_plot,
  './figures/fig6.png',
  units = 'in',
  width = 3.5,
  height = 6
)


# 6. CT data prep --------------------------------------------------------------

# Pivot wider so each row is a biological replicate
ct_values_wide = ct_values %>%
  mutate(
    id = paste(
      genotype,
      age,
      biological_replicate,
      sep = '_'
      )
    ) %>%
  mutate(id = gsub(' ', '_', id)) %>%
  pivot_wider(
    names_from = CT_replicate,
    names_prefix = 'CT_replicate_',
    values_from = CT_value
  )


# Calculate each BR's average CT value
row_means = rowMeans(
  ct_values_wide[,6:8],
  na.rm = T
  )

ct_values_wide_avg = ct_values_wide %>%
  mutate(row_means = row_means)

# Export
write.csv(
  ct_values_wide_avg,
  './data/ct_values_wide_avg.csv',
  row.names = F
  )


# In excel, for each BR, ΔCT=subtract the Actin 2 Average CT value from experimental gene Average CT value.
# Import excel calculations
delta_ct_wide = read.csv(
  './data/delta_ct_wide.csv',
  stringsAsFactors = T
  )


# 7. Compute ΔCT calibration values using WT as reference ----------------------

day11_calib_values = delta_ct_wide %>%
  filter(age == 'Day 11' & genotype == 'WT')

CRU3_day11_calib_mean = day11_calib_values %>%
  filter(gene == 'CRU3') %>%
  summarise(CRU3_day11_calib_mean = mean(delta_CT)) %>%
  unlist(use.names = F)

OLE1_day11_calib_mean = day11_calib_values %>%
  filter(gene == 'OLE1') %>%
  summarise(OLE1_day11_calib_mean = mean(delta_CT)) %>%
  unlist(use.names = F)


day13_calib_values = delta_ct_wide %>%
  filter(age == 'Day 13' & genotype == 'WT')

CRU3_day13_calib_mean = day13_calib_values %>%
  filter(gene == 'CRU3') %>%
  summarise(CRU3_day13_calib_mean = mean(delta_CT)) %>%
  unlist(use.names = F)

OLE1_day13_calib_mean = day13_calib_values %>%
  filter(gene == 'OLE1') %>%
  summarise(OLE1_day13_calib_mean = mean(delta_CT)) %>%
  unlist(use.names = F)


# Add calibration values back to df
delta_ct_calibrated = delta_ct_wide %>%
  mutate(calibration_value = case_when(
    age == 'Day 11' & gene == 'CRU3' ~ CRU3_day11_calib_mean,
    age == 'Day 13' & gene == 'CRU3' ~ CRU3_day13_calib_mean,
    age == 'Day 11' & gene == 'OLE1' ~ OLE1_day11_calib_mean,
    age == 'Day 13' & gene == 'OLE1' ~ OLE1_day13_calib_mean
    )
  )


# 8. Final ΔΔCT and fold change calculations -----------------------------------

final_calculations = delta_ct_calibrated %>%
  # ΔΔCT = (calibrator ΔCT) - (sample ΔCT)
  mutate(delta_delta_CT = delta_CT - calibration_value) %>%
  # fold change = 2^-ΔΔCT
  mutate(foldchange = 2^-delta_delta_CT)

# Export
write.csv(
  final_calculations,
  './data/sample_l2fc.csv',
  row.names = F
  )


# 9. Mean fold changes for each experimental group -----------------------------

# CRU3 Day 11
CRU3_WT_11 = final_calculations %>%
  filter(
    age == 'Day 11' &
    genotype == 'WT' &
    gene == 'CRU3'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )


CRU3_agp31_11 = final_calculations %>%
  filter(
    age == 'Day 11' &
    genotype == 'agp31 null' &
    gene == 'CRU3'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )

# CRU3 Day 13
CRU3_WT_13 = final_calculations %>%
  filter(
    age == 'Day 13' &
    genotype == 'WT' &
    gene == 'CRU3'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )

CRU3_agp31_13 = final_calculations %>%
  filter(
    age == 'Day 13' &
    genotype == 'agp31 null' &
    gene == 'CRU3'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )


# OLE1 Day 11
OLE1_WT_11 = final_calculations %>%
  filter(
    age == 'Day 11' &
    genotype == 'WT' &
    gene == 'OLE1'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )


OLE1_agp31_11 = final_calculations %>%
  filter(
    age == 'Day 11' &
    genotype == 'agp31 null' &
    gene == 'OLE1'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )


# OLE1 Day13
OLE1_WT_13 = final_calculations %>%
  filter(
    age == 'Day 13' &
    genotype == 'WT' &
    gene == 'OLE1'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )


OLE1_agp31_13 = final_calculations %>%
  filter(
    age == 'Day 13' &
    genotype == 'agp31 null' &
    gene == 'OLE1'
    ) %>%
  summarise(
    age = age,
    genotype = genotype,
    gene = gene,
    n = n(),
    mean_l2fc = mean(foldchange),
    SEM_l2fc = sd(foldchange)/sqrt(n)
    )


# row bind into one df
table_summary_statistics = rbind(
  CRU3_WT_11,
  CRU3_agp31_11,
  CRU3_WT_13,
  CRU3_agp31_13,
  OLE1_WT_11,
  OLE1_agp31_11,
  OLE1_WT_13,
  OLE1_agp31_13
  ) %>%
  distinct() %>%
  select(gene, genotype, age, mean_l2fc, SEM_l2fc, n) %>%
  arrange(gene, genotype) %>%
  mutate(
    mean_l2fc = round(mean_l2fc, 2),
    SEM_l2fc = round(SEM_l2fc, 2)
  )

write.csv(
  table_summary_statistics,
  './figures/table_summary_statistics.csv',
  row.names = F
  )


# 10. Plot mean fold changes in gene expression --------------------------------

expression_plot = ggplot(table_summary_statistics) +
  geom_col(
    aes(
      x = genotype,
      y = mean_l2fc,
      fill = gene
      ),
    alpha = 0.8
    ) +
  geom_errorbar(
    aes(
      ymin = mean_l2fc - SEM_l2fc,
      ymax = mean_l2fc + SEM_l2fc,
      x = genotype
      ),
    width = 0.2
    ) +
  scale_fill_manual(
    values = c('#C4AF9A', '#81523F')
    ) +
  facet_grid(gene ~ age) +
  labs(
    x = 'Genotype',
    y = 'Mean fold-change in expression'
    ) +
  theme(legend.position = 'none')

# Export
ggsave(
  expression_plot,
  './figures/fig7.png',
  units = 'in',
  width = 5,
  height = 4
)


# 11. Hypothesis testing done in Excel -----------------------------------------
