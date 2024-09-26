source('../universal_functions/FAXP_init.R')

# 1.SF1 --------------------------------------------------------------
my_proj <- c('ProteomEx (5.6 nL)', 'Current (1.55 nL)')

df1 <- rio::import('data/Supp Fig1.xlsx')
tblS1AB <- df1[c(1, 2, 3, 6)] %>%
  setNames(c('Project', 'Organism', '# peptides', '# proteins'))
tblS1AB[tblS1AB$Project == 'v1', 'Project'] <- my_proj[1]
tblS1AB[tblS1AB$Project == 'v2', 'Project'] <- my_proj[2]
tblS1AB$Project %<>% factor(levels = my_proj)
dfbar <- tblS1AB %>% pivot_longer(cols = -(Project:Organism), names_to = 'Type', values_to = 'Number')

plotS1AB <- ggplot(dfbar, aes(x = Project, y = Number, color = Project))+
  facet_wrap('Type', scales = 'free_y'
             ) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(my_proj),
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF1AB.pdf', plotS1AB, width = 6, height = 6)


data_anchor <- rio::import('data/PeptideYield_AnchorTime.xlsx')
tblS1C <- data_anchor
set.seed(10)
plotS1C <- ggplot(data_anchor, aes(x = `Anchor time`, y = `Peptide yield (μg)`, color = `Anchor time`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = 'Peptide yield (μg)')+
  scale_color_manual(values = c('#A07030', '#402060')) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'top') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(unique(data_anchor$`Anchor time`)),
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF1C.pdf', plotS1C, width = 6, height = 6)


# 2.SF3 --------------------------------------------------------------
#read data
info2 <- rio::import('data/Supp Fig3.xlsx')
pep_ls <- list.files('data/FASPv7_result_variable_DDA', '^peptide\\.tsv$', full.names = T, recursive = T)
pro_ls <- list.files('data/FASPv7_result_variable_DDA', '^protein\\.tsv$', full.names = T, recursive = T)
pep_ls %<>% str_subset(str_c(info2$`Name in the folder`, collapse = '|'))
pro_ls %<>% str_subset(str_c(info2$`Name in the folder`, collapse = '|'))
df2pep_ls <- lapply(pep_ls, rio::import)
df2pro_ls <- lapply(pro_ls, rio::import)

fileinfo2 <- data.frame(SampleName = c(pep_ls, pro_ls), Type = c(rep('pep', length(pep_ls)), rep('pro', length(pro_ls))))
fileinfo2$ID <- str_extract(fileinfo2$SampleName, str_c(info2$`Name in the folder`, '_\\d+', collapse = '|'))
fileinfo2$`Name in the folder` <- str_extract(fileinfo2$ID, str_c(info2$`Name in the folder`, collapse = '|'))
fileinfo2 %<>% left_join(info2)
names(df2pep_ls) <- fileinfo2$ID[fileinfo2$Type == 'pep']
names(df2pro_ls) <- fileinfo2$ID[fileinfo2$Type == 'pro']

#peptide matrix
df2pep <- plyr::ldply(df2pep_ls, .id = 'ID')
df2pep %<>% filter(!str_detect(Protein, '^CON_'))
df2pep %<>% filter(!str_detect(`Mapped Proteins`, '^CON_'))
mat2pep_count <- df2pep %>% pivot_wider(id_cols = ID, names_from = Peptide, values_from = `Spectral Count`) %>% column_to_rownames('ID') %>% t()

#protein matrix
df2pro <- plyr::ldply(df2pro_ls, .id = 'ID')
df2pro %<>% filter(!str_detect(Protein, '^CON_'))
df2pro %<>% filter(!str_detect(`Indistinguishable Proteins`, '^CON_'))
mat2pro_count <- df2pro %>% pivot_wider(id_cols = ID, names_from = `Protein ID`, values_from = `Total Spectral Count`) %>% column_to_rownames('ID') %>% t()


## 2.1 numbers (B) -----
# # peptides, # proteins
# % cysteine contained peptides, % carbamidomethyl mod
npep <- apply(mat2pep_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# peptides') %>% rownames_to_column('ID')
npro <- apply(mat2pro_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# proteins') %>% rownames_to_column('ID')
nC <- df2pep %>%
  group_by(ID) %>%
  summarise(`% cysteine contained peptides` = 100 * sum(str_detect(Peptide, 'C')) / length(Peptide))
nCmod <- df2pep %>% filter(str_detect(Peptide, 'C')) %>%
  group_by(ID) %>%
  summarise(`% carbamidomethyl modification` = 100 * sum(str_detect(`Assigned Modifications`, 'C')) / length(`Assigned Modifications`))


tblS3B <- fileinfo2 %>% select(-(SampleName:Type)) %>% distinct() %>% 
  left_join(npep) %>% left_join(npro) %>% left_join(nC) %>% left_join(nCmod)
tblS3B$`Name used in the image` %<>% factor(levels = info2$`Name used in the image`)
tblS3B %<>% arrange(`Name used in the image`)
dfbar <- tblS3B %>%
  pivot_longer(cols = -(ID:`Name used in the image`), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2), c(1,3), c(2,3)), function(indice){
  levels(tblS3B$`Name used in the image`)[indice]
})
plotS3B <- ggplot(dfbar, aes(x = `Name used in the image`, y = Number, color = `Name used in the image`))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF3B_pep_pro_C_Cmod.pdf', plotS3B, width = 10, height = 7)

## 2.2 % missed cleavage (C) -----
df_missCleav <- calc_missCleav(df2pep$Peptide)

tblS3C <- df2pep %>% left_join(df_missCleav %>% rename(Peptide = pepseq), relationship = 'many-to-many') %>%
  count(ID, MissedCleavage) %>%
  with_groups(ID, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  inner_join(tblS3B %>% select(ID:`Name used in the image`), .)

plotS3C <- ggplot(tblS3C, aes(x = `Name used in the image`, y = `ratio (%)`, color = `Name used in the image`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage)), scales = 'free_y',) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF')) +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF3C_missed_cleavage.pdf', plotS3C, width = 10, height = 7)


## 2.3 peptide length statistics (D) -----
tblS3D <- df2pep %>%
  inner_join(fileinfo2) %>% 
  distinct(`Name used in the image`, Peptide, `Peptide Length`) %>% 
  count(`Name used in the image`, `Peptide Length`)

plotS3D <- ggplot(tblS3D, aes(x = `Peptide Length`, y = n, color = `Name used in the image`, fill = `Name used in the image`)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = '', y = '', subtitle = 'Length of peptides')+
  # scale_y_continuous(labels = scales::scientific) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')
ggsave('SF3D_peptide_length.pdf', plotS3D, width = 5, height = 5)


# 3.SF4 --------------------------------------------------------------
info3 <- rio::import('data/Supp Fig4.xlsx')
pep_ls <- list.files('data/FASPv7_result_fixed_DDA', '^peptide\\.tsv$', full.names = T, recursive = T)
pro_ls <- list.files('data/FASPv7_result_fixed_DDA', '^protein\\.tsv$', full.names = T, recursive = T)
pep_ls %<>% str_subset(str_c(info3$`Name in the folder`, collapse = '|'))
pro_ls %<>% str_subset(str_c(info3$`Name in the folder`, collapse = '|'))
df3pep_ls <- lapply(pep_ls, rio::import)
df3pro_ls <- lapply(pro_ls, rio::import)

fileinfo3 <- data.frame(SampleName = c(pep_ls, pro_ls), Type = c(rep('pep', length(pep_ls)), rep('pro', length(pro_ls))))
fileinfo3$ID <- str_extract(fileinfo3$SampleName, str_c(info3$`Name in the folder`, '_\\d+', collapse = '|'))
fileinfo3$`Name in the folder` <- str_extract(fileinfo3$ID, str_c(info3$`Name in the folder`, collapse = '|'))
fileinfo3 %<>% left_join(info3)
names(df3pep_ls) <- fileinfo3$ID[fileinfo3$Type == 'pep']
names(df3pro_ls) <- fileinfo3$ID[fileinfo3$Type == 'pro']

#peptide matrix
df3pep <- plyr::ldply(df3pep_ls, .id = 'ID')
df3pep %<>% filter(!str_detect(Protein, '^CON_'))
df3pep %<>% filter(!str_detect(`Mapped Proteins`, '^CON_'))
mat3pep_count <- df3pep %>% pivot_wider(id_cols = ID, names_from = Peptide, values_from = `Spectral Count`) %>% column_to_rownames('ID') %>% t()

#protein matrix
df3pro <- plyr::ldply(df3pro_ls, .id = 'ID')
df3pro %<>% filter(!str_detect(Protein, '^CON_'))
df3pro %<>% filter(!str_detect(`Indistinguishable Proteins`, '^CON_'))
mat3pro_count <- df3pro %>% pivot_wider(id_cols = ID, names_from = `Protein ID`, values_from = `Total Spectral Count`) %>% column_to_rownames('ID') %>% t()


## 3.1 numbers (B) -----
# # peptides, # proteins
# % cysteine contained peptides, % carbamidomethyl mod
npep <- apply(mat3pep_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# peptides') %>% rownames_to_column('ID')
npro <- apply(mat3pro_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# proteins') %>% rownames_to_column('ID')
nC <- df3pep %>%
  group_by(ID) %>%
  summarise(`% cysteine contained peptides` = 100 * sum(str_detect(Peptide, 'C')) / length(Peptide))
nCmod <- df3pep %>% filter(str_detect(Peptide, 'C')) %>%
  group_by(ID) %>%
  summarise(`% carbamidomethyl modification` = 100 * sum(str_detect(`Assigned Modifications`, 'C')) / length(`Assigned Modifications`))


tblS4B <- fileinfo3 %>% select(-(SampleName:Type)) %>% distinct() %>% 
  left_join(npep) %>% left_join(npro) %>% left_join(nC) %>% left_join(nCmod)
tblS4B$`Name used in the image` %<>% factor(levels = info3$`Name used in the image`)
tblS4B %<>% arrange(`Name used in the image`)
dfbar <- tblS4B %>%
  pivot_longer(cols = -(ID:`Name used in the image`), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2), c(1,3), c(2,3)), function(indice){
  levels(tblS4B$`Name used in the image`)[indice]
})
plotS4B <- ggplot(dfbar, aes(x = `Name used in the image`, y = Number, color = `Name used in the image`))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF4B_pep_pro_C_Cmod.pdf', plotS4B, width = 10, height = 5)

## 3.2 % missed cleavage (C) -----
df_missCleav <- calc_missCleav(df3pep$Peptide)

tblS4C <- df3pep %>% left_join(df_missCleav %>% rename(Peptide = pepseq), relationship = 'many-to-many') %>%
  count(ID, MissedCleavage) %>%
  with_groups(ID, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  inner_join(tblS4B %>% select(ID:`Name used in the image`), .)

plotS4C <- ggplot(tblS4C, aes(x = `Name used in the image`, y = `ratio (%)`, color = `Name used in the image`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage)), scales = 'free_y',) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF')) +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF4C_missed_cleavage.pdf', plotS4C, width = 10, height = 7)


## 3.3 peptide length statistics (D) -----
tblS4D <- df3pep %>%
  inner_join(fileinfo3) %>% 
  distinct(`Name used in the image`, Peptide, `Peptide Length`) %>% 
  count(`Name used in the image`, `Peptide Length`)

plotS4D <- ggplot(tblS4D, aes(x = `Peptide Length`, y = n, color = `Name used in the image`, fill = `Name used in the image`)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = '', y = '', subtitle = 'Length of peptides')+
  # scale_y_continuous(labels = scales::scientific) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')
ggsave('SF4D_peptide_length.pdf', plotS4D, width = 5, height = 5)

# 4.SF9 --------------------------------------------------------------
## 4.1 read data ------
info4 <- rio::import('../Figure5/data/20230219_Batch_design_CRC_demo.xlsx', sheet = 'Sheet2') %>% select(1:7)
df4pg <- rio::import('../Figure5/data/whole_library/Whole_lib_report.pg_matrix.tsv')
df4pr <- rio::import('../Figure5/data/whole_library/Whole_lib_report.pr_matrix.tsv')

colnames(df4pg) %<>% str_remove('^.+\\\\')
colnames(df4pr) %<>% str_remove('^.+\\\\')

info4 <- data.frame(
  SampleName = colnames(df4pg)[-(1:5)],
  Batch = str_extract(colnames(df4pg)[-(1:5)], 'b\\d+_(\\d+|pool)')
) %>% full_join(info4)

mat4pg <- df4pg %>%
  column_to_rownames('Protein.Ids') %>% 
  select(-(Protein.Group:First.Protein.Description)) %>%
  log2()

mat4pep <- df4pr %>%
  group_by(Stripped.Sequence) %>% 
  select(-(Protein.Group:Precursor.Id)) %>%
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>%
  log2()

## 4.2 pool QC ------
### 4.2.1 identity (A)-----
df4_count <- data.frame(`# proteins` = apply(mat4pg, 2, function(y) sum(!is.na(y))),
                        `# peptides` = apply(mat4pep, 2, function(y) sum(!is.na(y))), check.names = F) %>% 
  rownames_to_column('SampleName') %>% 
  inner_join(info4, .)

tblS9A <- dfbar <- df4_count %>%
  filter(str_detect(Batch, 'pool')) %>% 
  select(SampleName, Batch, `# proteins`, `# peptides`) %>% 
  pivot_longer(cols = -(SampleName:Batch), names_to = 'Type', values_to = 'Number') %>% 
  mutate(Label = 'Pool')

plotS9A <- ggplot(dfbar, aes(x = Label, y = Number))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'none')
ggsave('SF9A_pep_pro_number.pdf', plotS9A, width = 4, height = 5)

### 4.2.2 CV (B) ---------
tblS9B_protein <- df4pg_poolcv <- mat4pg %>% t() %>% as.data.frame() %>%
  rownames_to_column('SampleName') %>%
  inner_join(info4, .) %>%
  filter(str_detect(Batch, 'pool')) %>%
  select(-(ID:Random), -SampleName) %>% 
  pivot_longer(cols = -Batch, names_to = 'Protein', values_to = 'Intensity', values_drop_na = T) %>% 
  group_by(Protein) %>% 
  mutate(Intensity = 2^Intensity) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean, Label = 'Pool')

plotS9B_pro <- ggplot(df4pg_poolcv, aes(x = Label, y = CV))+
  geom_violin(fill = '#AAAAAA', color = NA, width = 1)+
  # stat_boxplot(geom = 'errorbar', color = '#000000', width = 0.05)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Protein group')+
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

tblS9B_peptide <- df4pep_poolcv <- mat4pep %>% t() %>% as.data.frame() %>%
  rownames_to_column('SampleName') %>%
  inner_join(info4, .) %>%
  filter(str_detect(Batch, 'pool')) %>%
  select(-(ID:Random), -SampleName) %>% 
  pivot_longer(cols = -Batch, names_to = 'Peptide', values_to = 'Intensity', values_drop_na = T) %>% 
  group_by(Peptide) %>% 
  mutate(Intensity = 2^Intensity) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean, Label = 'Pool')

plotS9B_pep <- ggplot(df4pep_poolcv, aes(x = Label, y = CV))+
  geom_violin(fill = '#AAAAAA', color = NA, width = 1)+
  # stat_boxplot(geom = 'errorbar', color = '#000000', width = 0.05)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Peptide')+
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

plotS9B <- ggpubr::ggarrange(plotS9B_pep, plotS9B_pro, nrow = 1, ncol = 2)
ggsave('SF9B_pep_pro_CV.pdf', plotS9B, width = 5, height = 5)


## 4.3 missed cleavage (DDA data) (E) ----
ddapeps4 <- list.files('data/CRC_DDA', pattern = '^peptide\\.tsv$', recursive = T, full.names = T)
names(ddapeps4) <- str_extract(ddapeps4, 'exp_\\d+')

df4pepdda <- plyr::ldply(ddapeps4, rio::import, .id = 'ID')
pepList <- df4pepdda %>% count(ID, Peptide) %>% 
  plyr::dlply('ID', function(dfsub){
    dfsub$Peptide
  })
df_missCleav <- plyr::ldply(pepList, calc_missCleav, .id = 'SampleName')
tblS9D <- df_missCleav %>%
  count(ID, MissedCleavage) %>%
  with_groups(ID, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  mutate(Label = 'ProteomEx_v2_CRC_DDA')

plotS9D <- ggplot(tblS9D, aes(x = Label, y = `ratio (%)`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage)), scales = 'free_y',) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF'))
ggsave('SF9D_missed_cleavage.pdf', plotS9D, width = 8, height = 4)


## 4.4 quantified proteins (DIA data) (F) ----
df4pg_sample <- mat4pg %>% t() %>% as.data.frame() %>%
  rownames_to_column('SampleName') %>%
  inner_join(info4, .) %>%
  filter(!str_detect(Batch, 'pool'))

proList_patient <- plyr::dlply(df4pg_sample, 'Patient', function(dfsub){
  mat_tmp <- dfsub %>%
    column_to_rownames('Batch') %>%
    select(-(SampleName:Random)) %>% 
    t() %>% as.data.frame()
  mat_tmp <- mat_tmp[apply(mat_tmp, 1, function(x) any(!is.na(x))), ]
  return(rownames(mat_tmp))
})

proList_slide <- df4pg_sample %>% #of patient1
  filter(Patient == 'P1') %>% 
  plyr::dlply('Slide', function(dfsub){
    mat_tmp <- dfsub %>%
      column_to_rownames('Batch') %>%
      select(-(SampleName:Random)) %>% 
      t() %>% as.data.frame()
    mat_tmp <- mat_tmp[apply(mat_tmp, 1, function(x) any(!is.na(x))), ]
    return(rownames(mat_tmp))
  })
names(proList_slide) %<>% c('1' = 'P1S1', '2' = 'P1S2', '3' = 'P1S3')[.]

proList_region <- df4pg_sample %>% #of patient1 slide1
  filter(Patient == 'P1', Slide == 1) %>% 
  mutate(Region = ifelse(str_detect(Region, 'C'), 'C', Region)) %>% 
  mutate(Region = factor(Region, levels = c('N', 'L', 'H', 'C'))) %>% 
  plyr::dlply('Region', function(dfsub){
    mat_tmp <- dfsub %>%
      column_to_rownames('Batch') %>%
      select(-(SampleName:Random)) %>% 
      t() %>% as.data.frame()
    mat_tmp <- mat_tmp[apply(mat_tmp, 1, function(x) any(!is.na(x))), ]
    return(rownames(mat_tmp))
  })
names(proList_region) %<>% str_c('P1S1', .)

pdf('SF9E_upset1.pdf', width = 5, height = 6)
print(
  UpSetR::upset(
    UpSetR::fromList(proList_patient),
    sets = rev(names(proList_patient)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on', order.by = 'degree', # order.by = 'freq',
    mainbar.y.label = '# proteins',
    point.size = 3, line.size = 1, set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5), mb.ratio = c(0.6, 0.4)
  ))
print(
  UpSetR::upset(
    UpSetR::fromList(proList_slide),
    sets = rev(names(proList_slide)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on', order.by = 'degree', # order.by = 'freq',
    mainbar.y.label = '# proteins',
    point.size = 3, line.size = 1, set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5), mb.ratio = c(0.6, 0.4)
  ))
graphics.off()

pdf('SF9E_upset2.pdf', width = 10, height = 6)
print(
  UpSetR::upset(
    UpSetR::fromList(proList_region),
    sets = rev(names(proList_region)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on', order.by = 'degree', # order.by = 'freq',
    mainbar.y.label = '# proteins',
    point.size = 3, line.size = 1, set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5), mb.ratio = c(0.6, 0.4)
  ))
graphics.off()



tblS9E <- proList_patient %>%
  append(proList_slide) %>%
  append(proList_region) %>% 
  create_dataframe_from_list()

## 4.5 PVCA ------



# 5.SF11 --------------------------------------------------------------
## 5.1 Fig.S11A gel-making ------
df5Apg <- rio::import('data/homogenization_condition_DIA/Homogenization_report.pg_matrix.tsv')
df5Apr <- rio::import('data/homogenization_condition_DIA/Homogenization_report.pr_matrix.tsv')
colnames(df5Apg) %<>% str_remove('^.+\\\\')
colnames(df5Apr) %<>% str_remove('^.+\\\\')

info5A <- data.frame(SampleName = colnames(df5Apg)[-(1:5)])
info5A$Label <- str_extract(info5A$SampleName, 'Gel\\d+')

mat5Apg <- df5Apg %>%
  column_to_rownames('Protein.Ids') %>% 
  select(-(Protein.Group:First.Protein.Description)) %>%
  log2()

mat5Apep <- df5Apr %>%
  group_by(Stripped.Sequence) %>% 
  select(-(Protein.Group:Precursor.Id)) %>%
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>%
  log2()

tblS11A <- data.frame(`# proteins` = apply(mat5Apg, 2, function(y) sum(!is.na(y))),
                      `# peptides` = apply(mat5Apep, 2, function(y) sum(!is.na(y))), check.names = F) %>% 
  rownames_to_column('SampleName') %>% 
  inner_join(info5A, .)

dfbar <- tblS11A %>%
  pivot_longer(cols = -(SampleName:Label), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)), function(indice){
  unique(tblS11A$Label)[indice]
})
plotS11A <- ggplot(dfbar, aes(x = Label, y = Number, color = Label))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF11A_pep_pro.pdf', plotS11A, width = 7, height = 5)


## 5.2 Fig.S11B LC gradient ------
df5Bpg <- rio::import('data/LC_gradients_DIA/LC_gradients_report.pg_matrix.tsv')
df5Bpr <- rio::import('data/LC_gradients_DIA/LC_gradients_report.pr_matrix.tsv')
colnames(df5Bpg) %<>% str_remove('^.+\\\\')
colnames(df5Bpr) %<>% str_remove('^.+\\\\')

info5B <- rio::import('data/Supp Fig7.xlsx')
info5B <- data.frame(
  SampleName = colnames(df5Bpg)[-(1:5)],
  `Name in the folder` = str_extract(colnames(df5Bpg)[-(1:5)], str_c(info5B$`Name in the folder`, collapse = '|')),
  check.names = F
) %>% inner_join(info5B)

mat5Bpg <- df5Bpg %>%
  column_to_rownames('Protein.Ids') %>% 
  select(-(Protein.Group:First.Protein.Description)) %>%
  log2()

mat5Bpep <- df5Bpr %>%
  group_by(Stripped.Sequence) %>% 
  select(-(Protein.Group:Precursor.Id)) %>%
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>%
  log2()

tblS11B <- data.frame(`# proteins` = apply(mat5Bpg, 2, function(y) sum(!is.na(y))),
                      `# peptides` = apply(mat5Bpep, 2, function(y) sum(!is.na(y))), check.names = F) %>% 
  rownames_to_column('SampleName') %>% 
  inner_join(info5B, .)

# add FWHM.Scans
df5Bstat <- rio::import('data/LC_gradients_DIA/LC_gradients_report.stats.tsv')
df5Bstat$File.Name %<>% str_remove('^.+\\\\')
tblS11B <- df5Bstat %>% select(File.Name, FWHM.Scans) %>%
  rename(SampleName = File.Name) %>% 
  inner_join(tblS11B, .)

dfbar <- tblS11B %>%
  pivot_longer(cols = -(SampleName:`Name used in the image`), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2)), function(indice){
  unique(tblS11B$`Name used in the image`)[indice]
})
plotS11B <- ggplot(dfbar, aes(x = `Name used in the image`, y = Number, color = `Name used in the image`))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF11B_pep_pro_FWHM.pdf', plotS11B, width = 7.5, height = 5)

source('../universal_functions/FAXP_init.R')


# 6.S9F -----------------------------------------------------

## 6.1 read data ------
#protein matrix
dfpg <- rio::import('../Figure5/data/whole_library/Whole_lib_report.pg_matrix.tsv')
colnames(dfpg) %<>% str_remove('^.+\\\\')

#protein info
protinfo <- dfpg %>% select(Protein.Group:First.Protein.Description)

#protein matrix (log2 transform)
mat <- dfpg %>%
  select(-(Protein.Ids:First.Protein.Description)) %>%
  column_to_rownames('Protein.Group') %>% log2()

#sample info
info <- rio::import('../Figure5/data/20230219_Batch_design_CRC_demo.xlsx', sheet = 'Sheet2') %>% select(-ncol(.))
info <- data.frame(SampleName = colnames(mat),
                   Batch = str_extract(colnames(mat), 'b\\d+_(\\d+|pool)')) %>% 
  full_join(info)
info$BatchHead <- str_remove(info$Batch, '_.+$')
info$Slide %<>% factor()
info$Region <- ifelse(str_detect(info$Region, 'PC|CC'), 'C', info$Region)
info$Region %<>% factor(levels = c('N', 'L', 'H', 'C'), ordered = T)

#rename files in mat and regenerate dfpg
colnames(mat) %<>% str_extract('b\\d+_(\\d+|pool)')
dfpg <- mat %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>% 
  inner_join(info, .)



## 6.2 Data preparation ---------
protein_expression_matrix <- dfpg %>%
  select(-SampleName, -ID, -(Rep:BatchHead)) %>% 
  column_to_rownames('Batch') %>% 
  mutate(Slide = as.character(Slide)) %>% 
  filter(!is.na(Patient), !is.na(Slide))
protein_expression_matrix[which(protein_expression_matrix$Patient == 'P2'), 'Slide'] <- 4
protein_expression_matrix[which(protein_expression_matrix$Patient == 'P3'), 'Slide'] <- 5

log2_data <- protein_expression_matrix[, -(1:3)] %>% as.matrix()

missing <- protein_expression_matrix %>%
  group_by(Region) %>% 
  select(-Patient, -Slide) %>%
  summarise_all(function(y) sum(is.na(y)) / length(y)) %>%
  mutate(Region = as.character(Region)) %>% 
  rbind(c(Region = 'ALL', apply(log2_data, 2, function(x) sum(is.na(x)) / length(x))))

tmp <- missing %>%
  # filter(Region != 'ALL') %>% 
  column_to_rownames('Region') %>%
  # apply(2, function(y) sum(y < 0.2) >= 3 ) %>% 
  apply(2, function(y) any(y < 0.5))
table(tmp)
# FALSE  TRUE 
# 1084  5667 

log2_data_50na <- log2_data[, names(tmp[tmp])]
nafill <- min(log2_data_50na, na.rm = T) + log2(0.8)
log2_data_50na[is.na(log2_data_50na)] <- nafill
pem_50na <- cbind(protein_expression_matrix[, 1:3], log2_data_50na)



## 6.3 PVCA -----
theDataMatrix <- t(log2_data_50na)
expInfo <- protein_expression_matrix[, 1:3]
threshold <- 0.9
pvca <- my_pvcaBatchAssess(theDataMatrix, expInfo, threshold)

# Average proportion
df_pvca <- data.frame(t(pvca$dat)) %>% 
  cbind(pvca$label) %>% 
  setNames(c('RandomEffectWtAveProp', 'EffectName')) %>% 
  arrange(RandomEffectWtAveProp) %>% 
  mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
  arrange(desc(RandomEffectWtAveProp)) %>% 
  mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)
plot_pvca <- ggplot(df_pvca, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
  geom_bar(stat = 'identity', color = 'white') +
  coord_polar(theta = 'y', start = 0) +
  theme(legend.position = 'none') +
  geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
  scale_fill_manual(values = pca_colors) +
  theme_void() +
  xlim(0.5, 2.5)

# every PC
df_prop <- data.frame(pvca$matrixWtProp)
df_prop_percent <- data.frame(t(apply(df_prop, 1, function(x) x / sum(x))))
df_prop_percent$PC <- str_glue('PC{1:nrow(df_prop_percent)}: ({round(100 * pvca$eigenData$values[1:nrow(df_prop_percent)] / sum(pvca$eigenData$values[1:nrow(df_prop_percent)]), 2)} %)')
df_prop_percent$PC %<>% factor(., levels = .)

plot_pvca_individual <- plyr::dlply(df_prop_percent, 'PC', function(dfsub){
  tbl <- dfsub %>% select(-PC) %>% 
    setNames(., str_replace(names(.), '\\.', ':')) %>% 
    t() %>% data.frame() %>% 
    setNames('RandomEffectWtAveProp') %>% 
    rownames_to_column('EffectName') %>% 
    arrange(RandomEffectWtAveProp) %>% 
    mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
    arrange(desc(RandomEffectWtAveProp)) %>% 
    mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)
  
  set.seed(1000)
  ggplot(tbl, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
    geom_bar(stat = 'identity', color = 'white') +
    coord_polar(theta = 'y', start = 0) +
    geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
    scale_fill_manual(values = pca_colors) +
    labs(subtitle = dfsub$PC) +
    theme_void() +
    theme(legend.position = 'none') +
    xlim(0.5, 2.5)
})
length(plot_pvca_individual) # 100
legend_pies <- ggpubr::get_legend(plot_pvca)
p <- ggpubr::ggarrange(plotlist = plot_pvca_individual, nrow = 10, ncol = 10) %>%
  ggpubr::ggarrange(legend_pies, widths = c(0.95, 0.05))
plot(cumsum(pvca$eigenData$values / sum(pvca$eigenData$values)), xlab = 'PC_n', ylab = 'Summed eigenvalues proportion')


# top3 PCs
pvca$eigenData$values[1:3]
dfsub <- apply(df_prop_percent[1:3, 1:7], 2, function(y) {
  y * pvca$eigenData$values[1:3]
}) %>% as.data.frame() %>%
  mutate_all(sum) %>% slice(1)
dfsub <- dfsub / sum(dfsub)
tblS9F <- dfsub %>%
  setNames(., str_replace(names(.), '\\.', ':')) %>% 
  t() %>% data.frame() %>% 
  setNames('RandomEffectWtAveProp') %>% 
  rownames_to_column('EffectName') %>% 
  arrange(RandomEffectWtAveProp) %>% 
  mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
  arrange(desc(RandomEffectWtAveProp)) %>% 
  mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)


pS9F <- ggplot(tblS9F, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
  geom_bar(stat = 'identity', color = 'white') +
  coord_polar(theta = 'y', start = 0) +
  geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
  scale_fill_manual(values = pca_colors) +
  labs(subtitle = dfsub$PC) +
  theme_void() +
  xlim(0.5, 2.5)

ggsave('Figure_S9F_variance_top3PCs.pdf', pS9F, width = 5, height = 4)




# 7.S10B -----------------------------------------------------
load('../Figure5/Figure_5DE_6A-E_S9C_S10A.RData') # based on Figure5 data
dim(DEA) # 17574    16
DEP %>% pull(Protein) %>% unique() %>% length() # 335 DEPs in total

## disease progress associated proteins -----
my_comparisons <- c('L-N', 'H-L', 'C-H')
DEA_prog <- DEA %>% filter(Comparison %in% my_comparisons,
                           adj.P < 0.05)
prot_prog <- DEA_prog %>%
  pivot_wider(id_cols = 'Protein', names_from = 'Comparison', values_from = 'Log2FC') %>%
  column_to_rownames('Protein') %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x)) %>% # do not consider NA fold-change
  apply(1, function(x) (all(x >= 0) | all(x <= 0)) & any(abs(x) > log2(2))) %>% # fold-change values should be in the same trend, and at least one of them >2
  which() %>% names()
DEP_prog <- DEA_prog %>% filter(Protein %in% prot_prog)

prot_prog1 <- DEP_prog %>% count(Protein) %>% filter(n > 1) %>% pull(Protein) # at least two significant comparison
DEP_prog1 <- DEP_prog %>% filter(Protein %in% prot_prog1)

prot_prog2 <- DEP_prog1 %>% select(matches('^Missing.+_')) %>% apply(1, function(x) sum(x < 0.2) >= 3) %>% filter(DEP_prog1, .) %>% pull(Protein) %>% unique() # Missing <20% in at least three Regions
DEP_prog2 <- DEP_prog1 %>% filter(Protein %in% prot_prog2)


## heatmap ---------------------
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

clust_col <- hclust(dist(t(mat_heat[unique(DEP_prog2$Protein), ]), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()
heat3 <- pheatmap(mat_heat[unique(DEP_prog2$Protein), ],
                  annotation_col = dfheat %>% select(Region),
                  annotation_colors = anno_colors,
                  cluster_rows = T, cluster_cols = clust_col,
                  clustering_distance_rows = 'euclidean',
                  clustering_distance_cols = 'euclidean',
                  clustering_method = 'ward.D2',
                  show_rownames = F, show_colnames = F
)

bk <- unique(c(seq(-2, 2, length = 50)))
X <- mat_heat[heat3$tree_row$labels[heat3$tree_row$order], ]
rownames(X) <- DEP_prog2 %>% distinct(Protein, Label) %>% set_rownames(.$Protein) %>% .[rownames(X), ] %>% pull(Label)
pSF10B <- pheatmap(X,
                   annotation_col = dfheat %>% select(Region, Patient),
                   annotation_colors = anno_colors,
                   cluster_rows = F, cluster_cols = clust_col,
                   show_rownames = T, show_colnames = F,
                   color = heat_colors, breaks = bk,
                   fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                   border_color = F, na_col = '#AAAAAA',
                   filename = 'pSF10B_heatmap.pdf',
                   width = 10, height = 4
)
tblS10B <- X %>% as.data.frame() %>% rownames_to_column('Protein')


# Output ------------------------------------------------------------------
list(
  FigureS1AB = tblS1AB,
  FigureS1C = tblS1C,
  FigureS3B = tblS3B,
  FigureS3C = tblS3C,
  FigureS3D = tblS3D,
  FigureS4B = tblS4B,
  FigureS4C = tblS4C,
  FigureS4D = tblS4D,
  FigureS9A = tblS9A,
  FigureS9B_peptide = tblS9B_peptide,
  FigureS9B_protein = tblS9B_protein,
  FigureS9D = tblS9D,
  FigureS9E = tblS9E,
  FigureS9F = tblS9F,
  FigureS10B = tblS10B,
  FigureS11A = tblS11A,
  FigureS11B = tblS11B
) %>%
  rio::export('Source_data_SF1_3_4_9_10_11.xlsx')

save.image('Supplementary_figures.RData')



