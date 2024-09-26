source('../universal_functions/FAXP_init.R')


  
# 1.Read data and preprocess --------
#protein matrix
dfpro <- rio::import('data/reduction_alkylation_DDA/combined_protein.tsv')
dfpro %<>% filter(!str_detect(Protein, '^CON_'))
matpro <- dfpro %>%
  mutate(tmp = `Indistinguishable Proteins` %>%
           sapply(function(x){
             str_split(x, ', ')[[1]] %>%
               str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
               str_c(collapse = ', ')
           }),
         `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
  select(`Protein ID`, matches('^[^(MLQC)].+Total Spectral Count$')) %>%
  column_to_rownames('Protein ID')
colnames(matpro) %<>% str_remove(' Total Spectral Count$')

#sample info
info <- rio::import('data/F2_info.xlsx')
info <- data.frame(SampleName = colnames(matpro),
           'Label in file' = str_remove(colnames(matpro), '_\\d+$'),
           check.names = F) %>% left_join(info)

#peptide matrix
dfpep <- rio::import('data/reduction_alkylation_DDA/combined_peptide.tsv')
dfpep %<>% filter(!str_detect(Protein, '^CON_'))
matpep <- dfpep %>%
  select(`Peptide Sequence`, matches('^[^(MLQC)].+Spectral Count$')) %>%
  column_to_rownames('Peptide Sequence')
colnames(matpep) %<>% str_remove(' Spectral Count$')

#mod_peptide matrix
dfpep_mod <- rio::import('data/reduction_alkylation_DDA/combined_modified_peptide.tsv')
dfpep_mod %<>% filter(!str_detect(Protein, '^CON_'))
matpep_mod <- dfpep_mod %>%
  select(`Modified Sequence`, matches('^[^(MLQC)].+Spectral Count$')) %>%
  column_to_rownames('Modified Sequence')
colnames(matpep_mod) %<>% str_remove(' Spectral Count$')
# str_extract(dfpep_mod$`Modified Sequence`, '\\[[\\.\\d]+\\]') %>% unique()
# "[15.9949]" "[57.0214]" "[54.0474]" "[42.0106]"

# 2.Count numbers ---------------------------------------------------------
## 2.1 peptide number ---------
df <- matpep %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dflong <-  df %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Peptide', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA))
dfbar <- dflong %>% 
  group_by(SampleName) %>% 
  summarise(`# Peptides` = sum(!is.na(Count))) %>% 
  left_join(info)
lbls <- sort(unique(dfbar$`Label in figure`))
dfcomp <- data.frame(X = lbls) %>% cross_join(data.frame(X = lbls)) %>%
  filter(X.x < X.y)
dfcomp$p.value <- dfcomp %>% 
  t() %>% as.data.frame() %>% as.list() %>% 
  sapply(function(x){
    x1 <- dfbar %>% filter(`Label in figure` == x[1]) %>% pull(`# Peptides`)
    x2 <- dfbar %>% filter(`Label in figure` == x[2]) %>% pull(`# Peptides`)
    t.test(x1, x2)$p.value
  })
comparisons <- dfcomp %>% filter(p.value < 0.05) %>% select(-p.value) %>% t() %>% as.data.frame() %>% as.list()

set.seed(0)
p1 <- ggplot(dfbar, aes(x = `Label in figure`, y = `# Peptides`, color = `Label in figure`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '# peptides')+
  scale_color_brewer(palette = 'Set1') +
  scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,
    comparisons = comparisons,
    size = 5, hjust = 0.5, vjust = 0)

tbl1 <- dfbar %>% select(-Note) %>% as.data.frame() %>%
  cbind(NA) %>% cbind(NA) %>% 
  cbind(rbind(dfcomp, matrix(rep(NA, (nrow(dfbar) - nrow(dfcomp)) * ncol(dfcomp)), ncol = ncol(dfcomp)) %>% set_colnames(colnames(dfcomp))))
names(tbl1)[names(tbl1) == 'NA'] <- ''


## 2.2 protein number ---------
df <- matpro %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dflong <-  df %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Protein', values_to = 'Count') %>% 
  separate_rows('Protein') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA))
dfbar <- dflong %>% 
  group_by(SampleName) %>% 
  summarise(`# Proteins` = sum(!is.na(Count))) %>% 
  left_join(info)
lbls <- sort(unique(dfbar$`Label in figure`))
dfcomp <- data.frame(X = lbls) %>% cross_join(data.frame(X = lbls)) %>%
  filter(X.x < X.y)
dfcomp$p.value <- dfcomp %>% 
  t() %>% as.data.frame() %>% as.list() %>% 
  sapply(function(x){
    x1 <- dfbar %>% filter(`Label in figure` == x[1]) %>% pull(`# Proteins`)
    x2 <- dfbar %>% filter(`Label in figure` == x[2]) %>% pull(`# Proteins`)
    t.test(x1, x2)$p.value
  })
comparisons <- dfcomp %>% filter(p.value < 0.05) %>% select(-p.value) %>% t() %>% as.data.frame() %>% as.list()

set.seed(0)
p2 <- ggplot(dfbar, aes(x = `Label in figure`, y = `# Proteins`, color = `Label in figure`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '# proteins')+
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,
    comparisons = comparisons,
    size = 5, hjust = 0.5, vjust = 0)

tbl2 <- dfbar %>% select(-Note) %>% as.data.frame() %>%
  cbind(NA) %>% cbind(NA) %>% 
  cbind(rbind(dfcomp, matrix(rep(NA, (nrow(dfbar) - nrow(dfcomp)) * ncol(dfcomp)), ncol = ncol(dfcomp)) %>% set_colnames(colnames(dfcomp))))
names(tbl2)[names(tbl2) == 'NA'] <- ''



## 2.3 Mod carbamidomethyl ---------
df <- matpep_mod %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dflong <-  df %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'ModifiedPeptide', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA),
         Cmod57 = str_count(ModifiedPeptide, 'C\\[57\\.0214\\]'),
         Cwo57 = str_count(ModifiedPeptide, 'C[^\\[]'),
         Cmod57 = ifelse(is.na(Count), 0, Cmod57),
         Cwo57 = ifelse(is.na(Count), 0, Cwo57))

dfbar <- dflong %>% 
  group_by(SampleName) %>% 
  summarise(`% carbamidomethyl sites` = 100 * sum(Cmod57) / sum(Cmod57 + Cwo57)) %>% 
  left_join(info)
lbls <- sort(unique(dfbar$`Label in figure`))
dfcomp <- data.frame(X = lbls) %>% cross_join(data.frame(X = lbls)) %>%
  filter(X.x < X.y)
dfcomp$p.value <- dfcomp %>% 
  t() %>% as.data.frame() %>% as.list() %>% 
  sapply(function(x){
    x1 <- dfbar %>% filter(`Label in figure` == x[1]) %>% pull(`% carbamidomethyl sites`)
    x2 <- dfbar %>% filter(`Label in figure` == x[2]) %>% pull(`% carbamidomethyl sites`)
    t.test(x1, x2)$p.value
  })
comparisons <- dfcomp %>% filter(p.value < 0.05) %>% select(-p.value) %>% t() %>% as.data.frame() %>% as.list()

set.seed(0)
p3 <- ggplot(dfbar, aes(x = `Label in figure`, y = `% carbamidomethyl sites`, color = `Label in figure`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% carbamidomethyl sites')+
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,
    comparisons = comparisons,
    size = 5, hjust = 0.5, vjust = 0)

tbl3 <- dfbar %>% select(-Note) %>% as.data.frame() %>%
  cbind(NA) %>% cbind(NA) %>% 
  cbind(rbind(dfcomp, matrix(rep(NA, (nrow(dfbar) - nrow(dfcomp)) * ncol(dfcomp)), ncol = ncol(dfcomp)) %>% set_colnames(colnames(dfcomp))))
names(tbl3)[names(tbl3) == 'NA'] <- ''


## output -----------
X1 <- cbind(tbl1, tbl2, tbl3) %>% data.frame() %>%
  select(SampleName, Label.in.file, Label.in.figure,
         X..Peptides, X..Proteins, X..carbamidomethyl.sites) %>% 
  rename(`# Peptides` = X..Peptides,
         `# Proteins` = X..Proteins,
         `% carbamidomethyl sites` = X..carbamidomethyl.sites)

X2 <- cbind(tbl1, tbl2, tbl3) %>% data.frame() %>%
  select(matches('^X\\.[xy]$|^p\\.value')) %>% 
  rename(Group1 = X.x,
         Group2 = X.y,
         p.value.peptide = p.value,
         p.value.protein = p.value.1,
         p.value.carbamidomethyl.sites = p.value.2)
X <- X1 %>% cbind(NA) %>% cbind(NA) %>% cbind(X2)
names(X)[names(X) == 'NA'] <- ''
rio::export(X, 'F2I_source_data.xlsx')

p <- ggpubr::ggarrange(plotlist = list(p1, p2, p3), labels = c('I'), nrow = 1)
ggsave('F2I.pdf', p, width = 15, height = 5)


save.image('F2I.RData')
