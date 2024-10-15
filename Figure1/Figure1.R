source('../universal_functions/FAXP_init.R')
library(showtext)


# Read and process the data
df <- rio::import('data/timeline.xlsx')
df_proteomex <- df[1:10, 1:5] %>%
  set_colnames(.[1, ]) %>% slice(-1) %>% 
  mutate_at(vars(`T1 (hours)`:`Time lasting (hours)`), as.double) %>%
  rename(Process = ProteomEx)
df_faxp <- df[1:12, 7:11] %>%
  set_colnames(.[1, ]) %>% slice(-1) %>% 
  mutate_at(vars(`T1 (hours)`:`Time lasting (hours)`), as.double) %>%
  rename(Process = FAXP)

# Define progress levels
progress_all <- union(df_proteomex$Process, df_faxp$Process) %>%
  factor(levels = c(
    'Sectioning & de-waxing', 'Anchoring', 'Monomer incubation', 'Gelation',
    'Homogenization', 'Reduction & alkylation (whole gel)', 'Staining & expansion',
    'Imaging', 'Micro-dissection', 'Reduction and alkylation (punch), in-gel digestion & peptide extraction',
    'Filter-aided in-gel digestion & peptide elution', 'LC-MS/MS'
  )) %>%
  sort()

# Compute project-specific differences
proteomex_only <- setdiff(progress_all, df_faxp$Process)
faxp_only <- setdiff(progress_all, df_proteomex$Process)

# Prepare data for bar plots
dfbar <- list(ProteomEx = df_proteomex,
              FAXP = df_faxp) %>%
  plyr::ldply(.id = 'Project') %>%
  mutate(
    Process = fct_relevel(Process, rev(levels(progress_all))),
    Label = case_when(
      Process %in% proteomex_only ~ 'ProteomEx only',
      Process %in% faxp_only ~ 'FAXP only',
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(Project, desc(Process))

# Define numeric labels
num2chr <- setNames(c('①', '②', '③', '④', '⑤', '⑥', '⑦', '⑧', '⑨', '⑩', '⑪', '⑫'), 1:12)

# Prepare data for plotting
dfplot <- list(ProteomEx = df_proteomex, FAXP = df_faxp) %>%
  plyr::ldply(.id = 'Project') %>%
  mutate(
    Process = factor(Process, levels = levels(progress_all)),
    Label = case_when(
      Process %in% proteomex_only ~ 'ProteomEx only',
      Process %in% faxp_only ~ 'FAXP only',
      TRUE ~ NA_character_)
  ) %>%
  mutate(x = num2chr[as.numeric(Process)], .before = 1) %>% 
  group_by(x) %>%
  mutate(y2 = `Time lasting (hours)` / sum(`Time lasting (hours)`)) %>%
  ungroup() %>%
  arrange(Project, x)

# Fill missing values in `mat_procedure` using tidyr::fill
mat_procedure <- dfplot %>%
  pivot_wider(id_cols = Project, names_from = 'Process', values_from = `T2 (hours)`) %>% 
  column_to_rownames('Project') %>% 
  .[, levels(dfplot$Process)] %>% 
  add_column('0' = 0, .before = 1)
for(i in 1:nrow(mat_procedure)){
  for(j in 1:ncol(mat_procedure)){
    if(is.na(mat_procedure[i, j])){
      mat_procedure[i, j] <- mat_procedure[i, j-1]
    }
  }
}

dfplot_new <- mat_procedure %>% rownames_to_column('Project') %>% 
  pivot_longer(cols = -Project, names_to = 'Process', values_to = 'y1') %>%
  as.data.frame() %>% 
  mutate(
    Process = factor(Process, levels = c('0', as.character(progress_all))),
    x = ifelse(as.numeric(Process) != 1, num2chr[as.character(as.numeric(Process) - 1)], 0),
    x = factor(x, levels = c('0', num2chr))
  )

# Create the plots
p1 <- ggplot(dfplot_new, aes(x = x, y = y1, color = Project, fill = Project)) +
  geom_point(size = 2) +
  geom_line(aes(group = Project), linewidth = 1) +
  geom_text(data = dfplot_new %>% filter(y1 != 0),
            aes(label = y1), size = 5, vjust = -0.5, hjust = 0.5,
            show.legend = F) +
  labs(x = 'Procedures', y = 'Cumulative time (hours)') +
  scale_color_manual(values = project_colors) +
  theme_classic() +
  theme(text = element_text(size = 15, family = 'Arial Unicode MS'),
        axis.text.x.bottom = element_text(size = 20, family = 'Arial Unicode MS'),
        legend.position = 'top')
dfplot_lgd <- dfplot %>% distinct(x, Process) %>% arrange(x)
p2 <- ggplot(dfplot_lgd, aes(x = Process, color = Process, fill = x)) +
  geom_bar(width = 1) +
  scale_fill_manual(values = rep('white', nrow(dfplot_lgd))) +
  scale_color_manual(values = rep('white', nrow(dfplot_lgd))) +
  labs(x = '', y = 'Cumulative time (hours)') +
  theme_void() +
  theme(text = element_text(size = 15, family = 'Arial Unicode MS'),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'left')

p <- ggpubr::ggarrange(p1, p2)

# Save the plot
pdf('F1B_timeline.pdf', width = 13, height = 5)
showtext_begin()
print(p)
showtext_end()
graphics.off()

# Export the data and save the environment
rio::export(list(Figure1B = dfplot_new), 'F1B_timeline.xlsx')
save.image('F1B.RData')
