source('../universal_functions/FAXP_init.R')
library(showtext)


# read data -------------------------------------------------------------
df <- rio::import('data/timeline.xlsx')
df_proteomex <- df[1:10, 1:5] %>%
  set_colnames(.[1, ]) %>% slice(-1) %>% 
  mutate_at(vars(`T1 (hours)`:`Time lasting (hours)`), as.double)
df_faxp <- df[1:12, 7:11] %>%
  set_colnames(.[1, ]) %>% slice(-1) %>% 
  mutate_at(vars(`T1 (hours)`:`Time lasting (hours)`), as.double)
progress_all <- union(df_proteomex$ProteomEx, df_faxp$FAXP) %>% 
  factor(levels = c('Sectioning & de-waxing', 'Anchoring', 'Monomer incubation', 'Gelation', 'Homogenization', 'Reduction & alkylation (whole gel)', 'Staining & expansion', 'Imaging', 'Micro-dissection', 'Reduction and alkylation (punch), in-gel digestion & peptide extraction', 'Filter-aided in-gel digestion & peptide elution', 'LC-MS/MS')) %>% sort()
proteomex_only <- setdiff(progress_all, df_faxp$FAXP)
faxp_only <- setdiff(progress_all, df_proteomex$ProteomEx)


dfbar <- list(ProteomEx = df_proteomex %>% rename(Process = ProteomEx),
              FAXP = df_faxp %>% rename(Process = FAXP)) %>% 
  plyr::ldply(.id = 'Project') %>% 
  mutate(Process = factor(Process, levels = rev(levels(progress_all))))
dfbar[dfbar$Process %in% proteomex_only, 'Label'] <- 'ProteomEx only'
dfbar[dfbar$Process %in% faxp_only, 'Label'] <- 'FAXP only'
dfbar %<>% arrange(Project, desc(Process))


# number labels
num2chr <- c('①', '②', '③', '④', '⑤', '⑥', '⑦', '⑧', '⑨', '⑩', '⑪', '⑫') %>%
  setNames(1:12)
dfplot <- list(ProteomEx = df_proteomex %>% rename(Process = ProteomEx),
               FAXP = df_faxp %>% rename(Process = FAXP)) %>% 
  plyr::ldply(.id = 'Project') %>% 
  mutate(Process = factor(Process, levels = levels(progress_all))) %>% 
  add_column(x = as.numeric(.$Process) %>% num2chr[.], .before = 1)
dfplot[dfplot$Process %in% proteomex_only, 'Label'] <- 'ProteomEx only'
dfplot[dfplot$Process %in% faxp_only, 'Label'] <- 'FAXP only'
dfplot %<>% arrange(Project, x)

dfplot %<>% group_by(x) %>%
  reframe(Project = Project, y2 = `Time lasting (hours)` / sum(`Time lasting (hours)`)) %>%
  ungroup() %>% 
  left_join(dfplot, .)


# draw plot
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
  as.data.frame()
dfplot_new$Process %<>% factor(levels = c('0', as.character(progress_all)))
dfplot_new$x <- as.numeric(dfplot_new$Process) - 1
dfplot_new[dfplot_new$x != 0, 'x'] <- num2chr[dfplot_new[dfplot_new$x != 0, 'x']]
dfplot_new$x %<>% factor(levels = c('0', num2chr))

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
pdf('F1B_timeline.pdf', width = 13, height = 5)
showtext_begin()
print(p)
showtext_end()
graphics.off()

list(Figure1B = dfplot_new) %>% rio::export('F1B_timeline.xlsx')


save.image('F1B.RData')
