# Data analysis for Lab Report 1.


# 0. Load libraries & import data ----------------------------------------------
library(tidyverse)
library(ggpubr)

setwd("./Report1_SDSPAGE_WesternBlot")

# Bradford assay data
standards = read.csv(
  './data/standards_BSA_minus_empty.csv',
  header = T,
  stringsAsFactors = F
  )
samples = read.csv(
  './data/samples_minus_empty.csv',
  header = T,
  stringsAsFactors = F
  )

# CRU3 Molecular Weights
MW_STCurves_CRU3 = read.csv(
  './data/MW_STCurves_CRU3.csv',
  stringsAsFactors = F
  )

# HSP90 Molecular Weights
MW_STCurves_HSP90 = read.csv(
  './data/MW_STCurves_HSP90.csv',
  stringsAsFactors = F
  )


# 1. Data preparation ----------------------------------------------------------

# Pivot BSA standards and drop group 4 standards
standards_minus4 = standards[,-5]
standards_long = pivot_longer(standards_minus4, 2:4, names_to = 'group')


# Linear model of absorption as a function of concentration
std_lm = lm(data = standards_long, standards_long$value ~ standards_long$concentration_mg.mL)
summary(std_lm)
intercept = std_lm$coefficients[[1]]
slope = std_lm$coefficients[[2]]


# Find concentration of samples using lm
samples = samples[1:4, 1:2]
samples$concentration_mg.mL = (samples$value - intercept)/slope


# Combined data set with standards and samples
combined = full_join(samples, standards_long)
combined = combined %>%
  mutate(
    condition = case_when(
      !is.na(seed_type) ~ seed_type,
      is.na(seed_type) ~ 'BSA Standard',
    )
  ) %>%
  mutate(
    stock = concentration_mg.mL*2
  )

# Export
write.csv(combined, './figures/Table2_bradford_concentrations.csv', row.names = F)


# CRU3 Molecular Weights
licor_cru3 = MW_STCurves_CRU3 %>% filter(condition == 'Licor ladder')
cru3 = MW_STCurves_CRU3 %>% filter(condition != 'Licor ladder')


# HSP90 Molecular Weights
licor_HSP90 = MW_STCurves_HSP90 %>% filter(condition == 'Licor ladder')


# 2. Plot Bradford Assay Standard Curve ----------------------------------------

ggplot() +
  geom_smooth(
    aes(
      x = combined$concentration_mg.mL,
      y = combined$value
      ),
    method = 'lm',
    col = '#7f7f7f'
    ) +
  geom_point(
    aes(
      x = combined$concentration_mg.mL,
      y = combined$value,
      shape = combined$condition
      ),
    size = 3,
    alpha = 0.8
    ) +
  stat_regline_equation(
    aes(
      x = standards_long$concentration_mg.mL,
      y = standards_long$value,
      label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")
      )
    ) +
  scale_shape_manual(
    values = c(16, 1, 8, 15, 0),
    labels = c(
      'agp-31 null (Day 11)',
      'agp-31 null (Day 13)',
      'BSA Standard',
      'Ler-0 (Day 11)',
      'Ler-0 (Day 13)'
      )
    ) +
  labs(
    x = 'Protein Concentration (mg/mL)',
    y = 'Absorbance (at 595 nm)',
    title = 'Bradford Assay Standard Curve',
    shape = ''
    ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = 'right'
    )

# Export
ggsave('./figures/fig1.png', units = 'in', width = 7, height = 3.5)


# 3. Plot CRU3 Molecular Weight Standard Curve ---------------------------------

ggplot() +
  geom_smooth(
    aes(
      x = licor_cru3$px,
      y = licor_cru3$log_kD
      ),
    method = 'lm',
    col = '#7f7f7f'
    ) +
  geom_point(
    aes(
      x = MW_STCurves_CRU3$px,
      y = MW_STCurves_CRU3$log_kD,
      shape = MW_STCurves_CRU3$condition
      ),
    size = 3,
    alpha = 0.8
    ) +
  stat_regline_equation(
    aes(
      x = licor_cru3$px,
      y = licor_cru3$log_kD,
      label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")
      )
    ) +
  scale_shape_manual(
    values = c(1, 8),
    labels = c('CRU3-alpha', 'Licor Ladder')
    ) +
  labs(
    x = 'Band Size (px)',
    y = 'Molecular weight (log10(kD))',
    title = 'CRU3 Molecular Weight Standard Curve',
    shape = ''
    ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = 'bottom'
    )


# 4. Plot HSP90 Molecular Weight Standard Curve --------------------------------

ggplot() +
  geom_smooth(
    aes(
      x = licor_HSP90$px,
      y = licor_HSP90$log_kD
    ),
    method = 'lm',
    col = '#7f7f7f'
    ) +
  geom_point(
    aes(
      x = MW_STCurves_HSP90$px,
      y = MW_STCurves_HSP90$log_kD,
      shape = MW_STCurves_HSP90$condition
    ),
    size = 3,
    alpha = 0.8
    ) +
  stat_regline_equation(
    aes(
      x = licor_HSP90$px,
      y = licor_HSP90$log_kD,
      label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")
      )
    ) +
  scale_shape_manual(
    values = c(1, 8),
    labels = c('HSP90', 'Licor Ladder')
    ) +
  labs(
    x = 'Band Size (px)',
    y = 'Molecular weight (log10(kD))',
    title = 'HSP90 Molecular Weight Standard Curve',
    shape = ''
    ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = 'bottom'
    )
