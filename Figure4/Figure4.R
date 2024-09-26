source('../universal_functions/FAXP_init.R')


# Part I. Detection range -------
## 1.1 load data ----
#peptide matrix
dfpep1 <- rio::import('data/detection_range_DDA/combined_peptide.tsv')
dfpep1 %<>% filter(!str_detect(Protein, '^CON_'))
dfpep2 <- rio::import('data/detection_range_DIA/Detection_range_report.pr_matrix.tsv')

pep_dda <- dfpep1 %>%
  select(`Peptide Sequence`, matches('^In\\w+?_\\d+[mu]m_\\d+ Spectral Count$')) %>%
  column_to_rownames('Peptide Sequence')
colnames(pep_dda) %<>% str_remove(' Spectral Count$')

pep_dia <- dfpep2 %>%
  select(Stripped.Sequence, matches('InTip|InGel')) %>% 
  group_by(Stripped.Sequence) %>% 
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence')
colnames(pep_dia) %<>% str_remove('^.+\\\\')

#protein matrix
dfpro1 <- rio::import('data/detection_range_DDA/combined_protein.tsv')
dfpro1 %<>% filter(!str_detect(Protein, '^CON_'))
dfpro2 <- rio::import('data/detection_range_DIA/Detection_range_report.pg_matrix.tsv')

pro_dda <- dfpro1 %>%
  mutate(tmp = `Indistinguishable Proteins` %>%
           sapply(function(x){
             str_split(x, ', ')[[1]] %>%
               str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
               str_c(collapse = ', ')
           }),
         `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
  select(`Protein ID`, matches('^In\\w+?_\\d+[mu]m_\\d+ Total Spectral Count$')) %>%
  column_to_rownames('Protein ID')
colnames(pro_dda) %<>% str_remove(' Total Spectral Count$')

pro_dia <- dfpro2 %>%
  select(Protein.Group, matches('InTip|InGel')) %>% 
  column_to_rownames('Protein.Group')
colnames(pro_dia) %<>% str_remove('^.+\\\\')


#sample info
info <- rio::import('../Figure3/data/F3_info.xlsx', sheet = 2) %>%
  select(-Note) %>% 
  pivot_longer(cols = -(`Label in file (Diameter)`:`Label in file (Volume)`),
               values_to = 'Label in figure', values_drop_na = T) %>% 
  select(-name) %>%
  mutate(Method = str_remove(`Label in figure`, '_.+$'),
         Size = str_remove(`Label in figure`, '^.+_'))
info_dda <- data.frame(SampleName = colnames(pep_dda),
                       'Label in figure' = str_remove(colnames(pep_dda), '_\\d+$'),
                       MS = 'DDA',
                       check.names = F)
info_dia <- data.frame(SampleName = colnames(pep_dia),
                       'Label in figure' = str_extract(colnames(pep_dia), 'In\\w+?_\\d+[um]m'),
                       MS = 'DIA',
                       check.names = F)
info <- rbind(info_dda, info_dia) %>% left_join(info)
info %<>%
  mutate(Method = str_to_title(Method) %>% str_replace('In', 'In-'),
         MS = str_c(MS, '-MS'),
         Label = str_c(Method, ' ', MS))
rio::export(info, 'data/F4B_info.xlsx')

#regenerate dfpep and dfpro; within protein groups
dfpep_dda <- pep_dda %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpep_dia <- pep_dia %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro_dda <- pro_dda %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro_dia <- pro_dia %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)


## 1.2 Volume and peptide/protein counts (L) --------
dfline_pep <- c(apply(pep_dda, 2, function(y) sum(y != 0)), apply(pep_dia, 2, function(y) sum(!is.na(y)))) %>% 
  as.data.frame() %>% setNames('# peptides') %>% rownames_to_column('SampleName') %>% 
  inner_join(info, .)
dfline_pro <- c(apply(pro_dda, 2, function(y) sum(y != 0)), apply(pro_dia, 2, function(y) sum(!is.na(y)))) %>% 
  as.data.frame() %>% setNames('# proteins') %>% rownames_to_column('SampleName') %>% 
  inner_join(info, .)
tbl0 <- dfline <- dfline_pep %>% inner_join(dfline_pro)
# d_v_fct <- max(info$`Label in file (Diameter)`) / max(info$`Label in file (Volume)`)

dfline$`Label in file (Volume)` %<>% factor(levels = sort(unique(.)))
dfline$`Label in file (Diameter)` %<>% factor(levels = sort(unique(.)))
dfline_stat <- dfline %>%
  group_by(Label, Method, `Label in file (Volume)`, `Label in file (Diameter)`) %>%
  summarise_at(vars(`# peptides`, `# proteins`), list(mean = mean, sd = sd)) %>% 
  ungroup()


# guide functions
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}

# plot
plot_pep <- ggplot(dfline_stat)+
  aes(x = `Label in file (Volume)`, y = `# peptides_mean`, fill = Label) +
  geom_bar(stat = 'identity', position = position_dodge2(preserve = 'single'), width = 0.9, linewidth = 1)+
  geom_errorbar(aes(ymin = `# peptides_mean` - `# peptides_sd`, ymax = `# peptides_mean` + `# peptides_sd`), width = 0.9, position = position_dodge2(preserve = 'single'), linewidth = 0.5, color = '#000000') +
  geom_jitter(data = dfline, mapping = aes(y = `# peptides`), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75, seed = 15)) +
  labs(x = 'Volume (nL)', y = '# peptides', subtitle = '')+
  scale_fill_manual(values = detect_range_colors) +
  guides(x.sec = guide_axis_label_trans(~paste(as.character(sort(unique(info$`Label in file (Diameter)`)))))) +
  scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'right',
        axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x.top = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

plot_pro <- ggplot(dfline_stat)+
  aes(x = `Label in file (Volume)`, y = `# proteins_mean`, fill = Label) +
  geom_bar(stat = 'identity', position = position_dodge2(preserve = 'single'), width = 0.9, linewidth = 1)+
  geom_errorbar(aes(ymin = `# proteins_mean` - `# proteins_sd`, ymax = `# proteins_mean` + `# proteins_sd`), width = 0.9, position = position_dodge2(preserve = 'single'), linewidth = 0.5, color = '#000000') +
  geom_jitter(data = dfline, mapping = aes(y = `# proteins`), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75, seed = 15)) +
  labs(x = 'Volume (nL)', y = '# proteins', subtitle = '')+
  scale_fill_manual(values = detect_range_colors) +
  guides(x.sec = guide_axis_label_trans(~paste(as.character(sort(unique(info$`Label in file (Diameter)`)))))) +
  scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'right',
        axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x.top = element_text(angle = 0, vjust = 0.5, hjust = 0.5))


# Part I output -------
p4B <- ggpubr::ggarrange(plot_pep, plot_pro, nrow = 1, common.legend = T)
ggsave('F4B_detection_range.pdf', p4B, width = 11, height = 6)

list(Figure4B_detection_range = dfline_stat %>% mutate_at(vars(3:4), function(x) as.numeric(as.character(x))) %>% inner_join(tbl0, .)) %>%
  rio::export('F4_B_source_data.xlsx')


save.image('Figure4B.RData')


# Part II. Single cell/nucleus analysis --------
rm(list = ls())
gc()
source('../universal_functions/FAXP_init.R')

## 1.Read data and preprocess --------
# FAXP
df_faxp <- rio::import('data/Final_used_files_libfree_result_used/report.pg_matrix.tsv')
dim(df_faxp) # 4028 31
dfprot <- df_faxp %>% select(Protein.Group:First.Protein.Description)
df_faxp %<>%
  column_to_rownames('Protein.Group') %>% 
  select(-(Protein.Ids:First.Protein.Description)) %>% 
  setNames(str_remove(colnames(.), '^.+\\\\')) %>% 
  log2()
info <- data.frame(Sample = colnames(df_faxp))
info$Type <- str_extract(info$Sample, str_c(nucleus_type_order, collapse = '|'))
info$Type %<>% factor(., levels = nucleus_type_order)
info$ID <- str_extract(info$Sample, '\\d+\\.') %>% str_remove('\\.$') %>% str_c(info$Type, '_', .)
colnames(df_faxp) <- info$ID

info %<>% arrange(Type, ID)
df_faxp <- df_faxp[, info$ID]
comparisons <- combn(levels(info$Type)[1:3], 2, simplify = F)


# mouse_human_mapping
mouse_human_genes <- read.csv('data/informatics.jax.org/HOM_MouseHumanSequence.rpt',sep='\t')
mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
mouse <- mouse[,c(1,4)]
human <- human[,c(1,4)]
mh_data <- merge.data.frame(mouse, human, by = 'DB.Class.Key', all.y = T) # human list is longer than the mouse one
mh_data %<>% rename(Symbol.mouse = Symbol.x, Symbol.human = Symbol.y)


# human TFs
human_tf <- rio::import('data/humantfs/DatabaseExtract_v_1.01.xlsx')
human_tf %<>% filter(`Is TF?` == 'Yes')
dftf <- human_tf %>% inner_join(mh_data, by = c(`HGNC symbol` = 'Symbol.human'))

# MitoCarta3.0
mitocarta3 <- rio::import('data/mitocarta3/Mouse.MitoCarta3.0.xls', sheet = 3)

# blank
df_faxp_blk <- df_faxp %>% select(matches('blank')) %>%
  .[apply(., 1, function(x) any(!is.na(x))), ]
df_faxp_blk_top <- df_faxp_blk[apply(df_faxp_blk, 1, function(x) any(x[!is.na(x)] > log2(5e6))), ] # intensity > 5e6
dfprot %>%
  filter(Protein.Group %in% rownames(df_faxp_blk_top))



## 2.QC (F) ----------------------------------------------------------------
### 2.1 before QC filtering -----
# identity
prot_counter <- apply(df_faxp, 2, function(y) sum(!is.na(y)))
info$`# protein groups` <- prot_counter
info %<>% group_by(Type) %>%
  reframe(outlier_lower = get_outliers(`# protein groups`)[1],
          outlier_higher = get_outliers(`# protein groups`)[2]) %>% 
  inner_join(info, .)

set.seed(1)
pS7A <- ggplot(info, aes(x = Type, y = `# protein groups`, color = Type))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.3)+
  labs(x = '', y = '', subtitle = '# protein groups')+
  scale_color_manual(values = type_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
tblS7A <- info %>% arrange(Type, `# protein groups`)


# CV
df_cv <- 2^df_faxp %>% t() %>% as.data.frame() %>% rownames_to_column('ID') %>%
  inner_join(info %>% select(ID, Type), .) %>% 
  group_by(Type) %>% select(-ID) %>% summarise_all(cv)
pS7B <- df_cv %>% pivot_longer(cols = -Type, names_to = 'Protein.Group', values_to = 'CV') %>% 
  ggplot(aes(x = Type, y = CV))+
  geom_violin(fill = 'grey', color = NA) +
  geom_boxplot(width = 0.1) +
  labs(x = '', y = 'Coefficient of variation') +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

tblS7B <- df_cv %>%
  pivot_longer(cols = -Type, names_to = 'Protein.Group', values_to = 'CV')


# Pearson's correlation
mat_cor <- t(df_faxp)
rownames(mat_cor)

cors_sn <- df_faxp[, info$ID[info$Type == '1nucleus']] %>% 
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>% 
  round(2)
cors_c1n <- df_faxp[, info$ID[info$Type == 'cell_1nucleus']] %>% 
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>% 
  round(2)
cors_c2n <- df_faxp[, info$ID[info$Type == 'cell_2nucleus']] %>% 
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>% 
  round(2)

pdf('Figure_S7CDE.pdf', width = 8, height = 8)
my_corrplot(cors_sn)
my_corrplot(cors_c1n)
my_corrplot(cors_c2n)
graphics.off()

tblS7C <- cors_sn %>% as.data.frame() %>% rownames_to_column('ID')
tblS7D <- cors_c1n %>% as.data.frame() %>% rownames_to_column('ID')
tblS7E <- cors_c2n %>% as.data.frame() %>% rownames_to_column('ID')


### 2.2 QC filtering -----
bad_acquisitions <- c('1nucleus_3', 'cell_1nucleus_3', 'cell_1nucleus_1')
info %<>%
  filter(Type != 'blank_control', !(ID %in% bad_acquisitions)) %>% 
  mutate(Type = factor(Type, levels = levels(Type)[1:3]))
df_faxp <- df_faxp[, info$ID] %>% na_ratio_cutoff(cutoff = 1, MARGIN = 1)
dim(df_faxp) # 3998 18



### 2.3 after QC filtering (F) ------
# identity
set.seed(1)
p4F <- ggplot(info, aes(x = Type, y = `# protein groups`, color = Type))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.3)+
  labs(x = '', y = '# protein groups')+
  scale_color_manual(values = type_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none',
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  ggpubr::stat_compare_means(
    method = 't.test',
    comparisons = comparisons,
    size = 5, hjust = 0.5, vjust = 0)

tbl4F <- info %>% arrange(Type, `# protein groups`)


### 2.4 missing ratio ------
prot_miss <- apply(df_faxp, 1, function(x) sum(is.na(x)) / length(x))
sample_miss <- apply(df_faxp, 2, function(x) sum(is.na(x)) / length(x))
info$SampleMissingRatio <- sample_miss

# missing cutoff
df_miss_cutoff <- data.frame(
  `Missing ratio cutoff` = seq(0, 1, by = 0.05),
  `Remained protein number` = sapply(seq(0, 1, by = 0.05), function(cutoff){ sum(prot_miss <= cutoff) }),
  check.names = F
)
plot(df_miss_cutoff)

# cutoff == 50%
df_faxp05 <- df_faxp[prot_miss < 0.5, ]


### 2.5 Additional check ----------
#### location stat ---------
dfloc0 <- df_faxp %>%
  na_ratio_cutoff(cutoff = 1) %>% 
  rownames_to_column('Protein.Group') %>% 
  inner_join(dfprot, .)
dfloc0$Nucleus.Median <- dfloc0 %>%
  select(`1nucleus_1`:`1nucleus_6`) %>% 
  apply(1, function(x) log2(median(2^x, na.rm = T)))
dfloc0$CellMononuclear.Median <- dfloc0 %>%
  select(cell_1nucleus_2:cell_1nucleus_8) %>% 
  apply(1, function(x) log2(median(2^x, na.rm = T)))
dfloc0$CellBinuclear.Median <- dfloc0 %>%
  select(cell_2nucleus_1:cell_2nucleus_7) %>% 
  apply(1, function(x) log2(median(2^x, na.rm = T)))
dfloc0$All.Median <- dfloc0 %>%
  select(`1nucleus_1`:cell_2nucleus_7) %>% 
  apply(1, function(x) log2(median(2^x, na.rm = T)))

dfloc1 <- dfloc0 %>% 
  inner_join(., mitocarta3 %>% distinct(UniProt, MitoCarta3.0_Evidence, `HPA_Main_Location_2020 (Reliability)`, Tissues),
             by = c('Protein.Group' = 'UniProt'))
dfloc2 <- dfloc0 %>% anti_join(dfloc1) %>% 
  inner_join(., mitocarta3 %>% distinct(Symbol, MitoCarta3.0_Evidence, `HPA_Main_Location_2020 (Reliability)`, Tissues),
             by = c('Genes' = 'Symbol'))
dfloc3 <- dfloc0 %>% anti_join(dfloc1) %>% anti_join(dfloc2)
dfloc_whole <- rbind(dfloc1, dfloc2) %>% left_join(dfloc0, .) %>% 
  arrange(`HPA_Main_Location_2020 (Reliability)`)
x3 <- which(str_detect(dfloc_whole$`HPA_Main_Location_2020 (Reliability)`, ';'))
x2 <- which(is.na(dfloc_whole$`HPA_Main_Location_2020 (Reliability)`))
x1 <- 1:nrow(dfloc_whole) %>% setdiff(c(x2, x3))
dfloc_whole <- dfloc_whole[c(x1, x2, x3), ]

# dfloc <- dfloc_whole # NA as Others
dfloc <- rbind(dfloc1, dfloc2) %>%
  filter(!is.na(`HPA_Main_Location_2020 (Reliability)`)) %>% # omit NA
  filter(!str_detect(`HPA_Main_Location_2020 (Reliability)`, ';')) # remove multiple matching

# set "Others"
choice <- c('duplicate.split', 'as.others')[2]
proportion_min_cutoff <- c(0.01, 0.01)[2]
proportion_min_cutoff_lbl <- str_glue("<{proportion_min_cutoff*100}%")

if(choice == 'duplicate.split'){
  X_loc <- dfloc %>%
    rename(Location = `HPA_Main_Location_2020 (Reliability)`) %>% 
    select(matches('Median$'), Location, Protein.Group) %>% 
    mutate(Location = str_remove(Location, ' \\(.+\\)$')) %>% 
    mutate(Location = ifelse(is.na(Location), 'Others', Location)) %>% 
    # mutate(Location = ifelse(str_detect(Location, ';'), 'Others', Location)) %>%
    separate_rows(Location, sep = ';') %>% 
    pivot_longer(1:4, names_to = 'Type', values_to = 'MedianIntensity') %>% 
    mutate(Type = str_remove(Type, '\\.Median'),
           MedianIntensity = 2^MedianIntensity)
} else if (choice == 'as.others'){
  X_loc <- dfloc %>%
    rename(Location = `HPA_Main_Location_2020 (Reliability)`) %>% 
    select(matches('Median$'), Location, Protein.Group) %>% 
    mutate(Location = str_remove(Location, ' \\(.+\\)$')) %>% 
    mutate(Location = ifelse(is.na(Location), 'Others', Location)) %>% 
    mutate(Location = ifelse(str_detect(Location, ';'), 'Others', Location)) %>%
    pivot_longer(1:4, names_to = 'Type', values_to = 'MedianIntensity') %>% 
    mutate(Type = str_remove(Type, '\\.Median'),
           MedianIntensity = 2^MedianIntensity)
}

# stat
X_loc %>% count(Location) %>% arrange(desc(n))
stat_loc <- plyr::ddply(X_loc, 'Type', function(dfsub){
  dfsub %<>%
    select(Type, Location, MedianIntensity) %>% 
    group_by(Location) %>%
    summarise(SumMedInten = sum(MedianIntensity, na.rm = T))
  dfsub$Proportion <- dfsub$SumMedInten / sum(dfsub$SumMedInten, na.rm = T)
  return(dfsub)
})

# visual
dfbar_loc <- stat_loc %>% distinct(Type, Location, Proportion)
location_order <- dfbar_loc %>% filter(Type == 'All') %>% 
  arrange(desc(Proportion)) %>% 
  pull(Location) %>% unique() %>%
  rev() %>% c(proportion_min_cutoff_lbl, .) %>% unique() %>% rev()
locationAnno_order <- dfbar_loc %>%
  filter(Type == 'All') %>%
  mutate(Location = factor(Location, levels = location_order)) %>% 
  arrange(Location) %>% 
  mutate(LocationAnno = ifelse(Location == proportion_min_cutoff_lbl, Location, str_glue("{Location} ({round(Proportion*100, 1)} %)"))) %>% 
  mutate(LocationAnno = factor(LocationAnno, levels = unique(LocationAnno))) %>% 
  pull(LocationAnno) %>% as.character() %>% c(., proportion_min_cutoff_lbl)
locAnnoConvertor <- locationAnno_order %>% setNames(location_order)
dfbar_loc$Location[dfbar_loc$Proportion < proportion_min_cutoff] <- proportion_min_cutoff_lbl
dfbar_loc$Location <- locAnnoConvertor[dfbar_loc$Location]

dfbar_loc %<>%
  filter(Location == proportion_min_cutoff_lbl) %>% 
  group_by(Type, Location) %>% 
  summarise(Proportion = sum(Proportion)) %>% 
  as.data.frame() %>% 
  rbind(dfbar_loc %>% filter(Location != proportion_min_cutoff_lbl), .)
dfbar_loc %<>% mutate(
  Type = factor(Type, levels = c('Nucleus', 'CellMononuclear', 'CellBinuclear', 'All')),
  Location = factor(Location, levels = locationAnno_order)
) %>% arrange(Type, Location)
nLocation <- length(unique(dfbar_loc$Location))

set.seed(1)
pbar_location <- ggplot(dfbar_loc) +
  aes(x = Type, y = Proportion, fill = Location) +
  geom_col() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, 'Spectral'))(nLocation)[c(sample(seq(1, nLocation, 2)), sample(seq(2, nLocation, 2)))]) +
  labs(x = '') +
  ggthemes::theme_par()
ggsave('location_mitocarta3Anno.pdf', pbar_location, width = 10, height = 9)

list(`barplot source data` = dfbar_loc,
     `FAXP data-MitoCarta3 annotated` = dfloc_whole,
     `FAXP data stat` = X_loc %>% inner_join(stat_loc) %>%
       inner_join(dfprot, .) %>% rename(`LocSum-MedInten` = SumMedInten)
) %>% 
  rio::export('location_mitocarta3Anno.xlsx')



#### median protein intensity -----
## marking TFs and histones ##
mouse_tfs <- dfprot %>% filter(Genes %in% dftf$Symbol.mouse)
faxp_median_intensity <- function(faxp_matrix, pdfNameLabel = ''){
  # median intensity and labeling
  median_inten <- faxp_matrix %>% t() %>% as.data.frame() %>% 
    rownames_to_column('ID') %>%
    inner_join(info %>% select(ID, Type), .) %>% 
    group_by(Type) %>% 
    summarise_if(is.double, median, na.rm = T) %>% 
    rbind(faxp_matrix %>% t() %>% as.data.frame() %>% 
            summarise_all(median, na.rm = T) %>% 
            add_column(Type = 'All', .before = 1)) %>% 
    column_to_rownames('Type') %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('Protein.Group') %>% 
    inner_join(dfprot, .)
  median_inten$IsTF <- sapply(median_inten$Protein.Group, function(x) any(unlist(str_split(x, ';')) %in% mouse_tfs$Protein.Group))
  median_inten$IsHistone <- sapply(median_inten$Protein.Group, function(x) any(unlist(str_split(x, ';')) %in% mouse_histone$Entry))
  
  median_inten_all <- median_inten %>% 
    arrange(desc(All)) %>% 
    mutate(Rank = 1:nrow(.))
  
  # subsets
  histone_median <- median_inten_all %>%
    filter(IsHistone)
  tf_median <- median_inten_all %>%
    filter(IsTF)
  
  # visualization
  p_inten0 <- ggplot()+
    geom_point(data = subset(median_inten_all, !IsTF),
               aes(x = Rank, y = All), size = 2) +
    labs(x = 'Rank', y = 'Log2-MedianIntensity')+
    theme_classic() +
    theme(text = element_text(size = 15))
  p_lbl_tf <- 
    geom_text_repel( # TFs
      data = tf_median,
      aes(x = Rank, y = All, label = Genes), seed = 1,
      # direction = 'y',
      nudge_y = 0.5,
      point.padding = 0.5, max.overlaps = 100,
      size = 3, color = '#C81366',# fill = NA, # alpha = 0.9,
      segment.size = 0.5, segment.color = '#C81366',
    )
  p_lbl_histone <- 
    geom_text_repel( # histones
      data = histone_median,
      aes(x = Rank, y = All, label = str_remove(Genes, ';.+$')), seed = 10,
      # direction = 'x', nudge_x = 2,
      point.padding = 0.5, max.overlaps = 100,
      size = 3, color = '#2953E5',# fill = NA, alpha = 0.9,
      segment.size = 0.5, segment.color = '#2953E5',
    )
  
  
  p_inten_tf <- p_inten0 +
    p_lbl_tf +
    geom_point(data = subset(median_inten_all, IsTF),
               aes(x = Rank, y = All), color = '#C81366', size = 4)
  p_inten_histone <- p_inten0 +
    p_lbl_histone +
    geom_point(data = subset(median_inten_all, IsHistone),
               aes(x = Rank, y = All), color = '#2953E5', size = 4)
  p_inten <- p_inten0 +
    p_lbl_tf +
    p_lbl_histone +
    geom_point(data = subset(median_inten_all, IsTF),
               aes(x = Rank, y = All), color = '#C81366', size = 4) +
    geom_point(data = subset(median_inten_all, IsHistone),
               aes(x = Rank, y = All), color = '#2953E5', size = 4)
  
  ggsave(str_glue('FAXP_protein_median_intensity_split_{pdfNameLabel}.pdf'),
         ggpubr::ggarrange(p_inten_tf, p_inten_histone, nrow = 1),
         width = 20, height = 10)
  ggsave(str_glue('FAXP_protein_median_intensity_{pdfNameLabel}.pdf'),
         p_inten, width = 10, height = 10)
  rio::export(median_inten_all, str_glue('FAXP_protein_median_intensity_{pdfNameLabel}.xlsx'))
  
}
faxp_median_intensity(df_faxp05, pdfNameLabel = 'miss05')



## 3.nucleus / cell (Mono) -------
### 3.1 data preparation -----
infoX <- info %>%
  filter(Type %in% c('1nucleus', 'cell_1nucleus')) %>% 
  select(ID, Type)

prot_cell_spec <- rownames(df_faxp)[apply(df_faxp[, info$ID[info$Type == '1nucleus']], 1, function(x) all(is.na(x))) & apply(df_faxp[, info$ID[info$Type != '1nucleus']], 1, function(x) sum(is.na(x)) / length(x) < 0.5)]
X_cell_spec <- df_faxp[prot_cell_spec, ]

X_cell_spec_imp <- X_cell_spec
X_cell_spec_imp[, infoX$ID[infoX$Type == '1nucleus']] %<>%
  impute_min(minvalue = apply(X_cell_spec, 1, min, na.rm = T),
             downfactor = 0.9,
             SD = apply(X_cell_spec, 1, sd, na.rm = T), width = 0.3)
X <- df_faxp[!(rownames(df_faxp) %in% prot_cell_spec), infoX$ID] %>%
  na_ratio_cutoff(cutoff = 0.5) %>% 
  rbind(X_cell_spec_imp[, infoX$ID])

list(matrix = X %>% rownames_to_column('rownames'),
     info = infoX) %>%
  rio::export('data/nucleus_cell_data.xlsx')

### 3.2 Pearson (G) ---------
cors <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
pdf('F4G.pdf', width = 8, height = 8)
my_corrplot(cors)
graphics.off()
tbl4G <- cors %>% as.data.frame() %>% rownames_to_column('ID')


### 3.3 DEA ----------
dfX <- X %>% t() %>% as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  inner_join(infoX, .) %>% 
  select(-ID) %>% 
  pivot_longer(cols = -Type,values_drop_na = T,
               names_to = 'Protein', values_to = 'Log2Intensity') %>% 
  as.data.frame()
dea <- plyr::ddply(dfX, 'Protein', function(dfsub){
  type1 <- levels(dfsub$Type)[1]
  type2 <- levels(dfsub$Type)[2]
  x1 <- dfsub$Log2Intensity[dfsub$Type == type1]
  x2 <- dfsub$Log2Intensity[dfsub$Type == type2]
  fc <- mean(2^x1) / mean(2^x2)
  pvalue <- tryCatch(t.test(x1, x2, var.equal = F)$p.value,
                     error = function(e) NA)
  data.frame(Log2FC = log2(fc), P.Value = pvalue)
})
dea %<>% inner_join(dfprot %>% rename(Protein = Protein.Group), .)
dep <- dea %>%
  filter(abs(Log2FC) > log2(2)) %>%
  mutate(adj.P = p.adjust(P.Value, 'BH')) %>% 
  filter(adj.P < 0.05)
dim(dep) # 416 8


X_dep <- X[dep$Protein, ]
pheatmap(X_dep, scale = 'row',
         color = heat_colors, na_col = '#AAAAAA',
         cluster_rows = T, cluster_cols = T,
         clustering_method = 'ward.D2',
         show_rownames = F, show_colnames = T
)



### 3.4 GSEA (H) ------------
enrichment_analysis_env_load()

dep_gse <- dea %>% filter(abs(Log2FC) > log2(1.5)) %>%
  mutate(adj.P = p.adjust(P.Value, 'BH')) %>% 
  filter(adj.P < 0.05)
tmp1 <- bitr(dep_gse$Genes, OrgDb = org.Mm.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID')
tmp2 <- bitr(dep_gse$Protein, OrgDb = org.Mm.eg.db, fromType = 'UNIPROT', toType = 'ENTREZID')
tmp <- tmp1 %>% full_join(tmp2) %>% setNames(c('Genes', 'ENTREZID', 'Protein'))
tmp %>% count(Genes) %>% drop_na() %>% count(n) # each n equals 1, except Tmpo 
dep_gse %<>% left_join(tmp) %>% arrange(desc(Log2FC))


gene_list <- drop_na(dep_gse, Genes)$Log2FC %>% setNames(drop_na(dep_gse, Genes)$Genes) %>% .[!duplicated(names(.))]
gid_list <- drop_na(dep_gse, ENTREZID)$Log2FC %>% setNames(drop_na(dep_gse, ENTREZID)$ENTREZID) %>% .[!duplicated(names(.))]


#### GOCC -----
gseGOCC <- gseGO(gene_list, OrgDb = org.Mm.eg.db,
                 ont = 'CC', keyType = 'SYMBOL',
                 minGSSize = 10, maxGSSize = 300,
                 exponent = 1, eps = 0,
                 pvalueCutoff = 1, pAdjustMethod = 'BH',
                 verbose = T, seed = F, by = 'fgsea'
)
tbl4H1 <- gseGOCC[]
interseted_cc <- c('nucleoplasm', 'mitochondrion', 'spliceosomal complex', 'cytoplasm', 'nucleolus', 'protein-DNA complex', 'ribonucleoprotein complex', 'chromosome', 'endoplasmic reticulum', 'Golgi apparatus', 'peroxisome')

p4H1 <- gseGOCC %>%
  filter(Description %in% interseted_cc) %>% 
  #filter(p.adjust < 5e-2, abs(NES) > 1.8) %>%
  ridgeplot(showCategory = 20, orderBy = 'pvalue', decreasing = T)+
  scale_fill_viridis_c(option = 'B', trans = 'log10',
                       breaks = c(1e-13, 1e-10, 1e-7, 1e-4, 5e-2),
                       limits = c(1e-13, 5e-2)) +
  labs(x = 'Log2FoldChange', subtitle = 'GSEA (GOCC)') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))


#### Reactome -------
gseReactome <- gsePathway(gid_list, organism = 'mouse',
                          minGSSize = 5, maxGSSize = 300,
                          exponent = 1, eps = 0,
                          pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                          verbose = T)
tbl4H2 <- gseReactome[]
interseted_rct <- c("Metabolism of RNA", "Processing of Capped Intron-Containing Pre-mRNA", "mRNA Splicing - Major Pathway", "Major pathway of rRNA processing in the nucleolus and cytosol", "Gene expression (Transcription)", "mRNA 3'-end processing", "Metabolism", "RNA Polymerase II Transcription", "Processing of Capped Intronless Pre-mRNA", "The citric acid (TCA) cycle and respiratory electron transport", "Transport of Mature Transcript to Cytoplasm", "DNA Double-Strand Break Repair")
p4H2 <- gseReactome %>%
  filter(p.adjust < 5e-2, abs(NES) > 1.8) %>%
  filter(Description %in% interseted_rct) %>% 
  ridgeplot(showCategory = 20, orderBy = 'pvalue', decreasing = T)+
  scale_fill_viridis_c(option = 'B', trans = 'log10', breaks = c(1e-13, 1e-10, 1e-7, 1e-4, 5e-2), limits = c(1e-13, 5e-2)) +
  labs(x = 'Log2FoldChange', subtitle = 'GSEA (Reactome)') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))



## 4.cell (Mono) / cell (Bi) -------
source('../universal_functions/FAXP_init.R')


### 4.1 Profile -----------
infoX <- info %>%
  filter(Type %in% c('cell_1nucleus', 'cell_2nucleus')) %>% 
  select(ID, Type)
X0 <- df_faxp[, infoX$ID] %>% na_ratio_cutoff(cutoff = 1)


# upset
dfX0 <- X0 %>% t() %>% as.data.frame() %>% 
  rownames_to_column('ID') %>% inner_join(infoX, .) %>% 
  select(-ID) %>% 
  pivot_longer(cols = -Type, values_drop_na = T,
               names_to = 'Protein', values_to = 'Log2Intensity') %>% 
  as.data.frame()
tmp <- dfX0 %>%
  distinct(Type, Protein) %>%
  mutate(Type = as.character(Type))
protein_detect_split <- split(tmp$Protein, tmp$Type) %>%
  setNames(c('Mono', 'Bi'))

graphics.off()
pdf('Figure_S8C.pdf', width = 10, height = 6)
print(
  UpSetR::upset(
    UpSetR::fromList(protein_detect_split),
    sets = names(protein_detect_split),
    keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on',
    order.by = 'freq',# order.by = c('freq', 'degree'),
    mainbar.y.label = '# proteins',
    point.size = 3, 
    line.size = 1,
    set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5),
    mb.ratio = c(0.6, 0.4)
  )
)
graphics.off()

protein_detect_venn <- list(`Mono-specific` = setdiff(protein_detect_split$Mono, protein_detect_split$Bi),
                            `Bi-specific` = setdiff(protein_detect_split$Bi, protein_detect_split$Mono),
                            overlap = intersect(protein_detect_split$Mono, protein_detect_split$Bi))
tblS8C <- create_dataframe_from_list(protein_detect_venn)


### 4.2 DEA -----------
X <- df_faxp[, infoX$ID] %>% na_ratio_cutoff(cutoff = 0.5)

dfX <- X %>% t() %>% as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  inner_join(infoX, .) %>% 
  select(-ID) %>% 
  pivot_longer(cols = -Type,values_drop_na = T,
               names_to = 'Protein', values_to = 'Log2Intensity') %>% 
  as.data.frame()
dea <- plyr::ddply(dfX, 'Protein', function(dfsub){
  type1 <- levels(dfsub$Type)[2]
  type2 <- levels(dfsub$Type)[3]
  x1 <- dfsub$Log2Intensity[dfsub$Type == type1]
  x2 <- dfsub$Log2Intensity[dfsub$Type == type2]
  fc <- mean(2^x1) / mean(2^x2)
  pvalue <- tryCatch(t.test(x1, x2, var.equal = F)$p.value,
                     error = function(e) NA)
  data.frame(Log2FC = log2(fc), P.Value = pvalue)
})
dea %<>% inner_join(dfprot %>% rename(Protein = Protein.Group), .)
dep <- dea %>%
  filter(abs(Log2FC) > log2(1.5)) %>%
  mutate(adj.P = p.adjust(P.Value, 'BH')) %>% 
  filter(adj.P < 0.05)
dim(dep) # 0  8

list(matrix = X %>% rownames_to_column('rownames'),
     info = infoX) %>%
  rio::export('data/cell_mono_bi_data.xlsx')

### 4.3 pearson's correlation ----
cors <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
pdf('Figure_S8A.pdf', width = 8, height = 8)
my_corrplot(cors)
graphics.off()

tblS8A <- cors %>% as.data.frame() %>% rownames_to_column('ID')



## 5.Histone of all ------
# read data
mouse_histone <- rio::import('data/uniprotkb_protein_name_Histone_AND_taxo_2024_08_07.xlsx')
mouse_histone %<>% filter(str_detect(`Protein names`, '^Histone '))
mouse_histone %<>% filter(str_detect(`Entry Name`, '^H\\d|^CENPA_MOUSE$'))

X1 <- rio::import('data/nucleus_cell_data.xlsx', sheet = 'matrix') %>% 
  column_to_rownames('rownames')
info1 <- rio::import('data/nucleus_cell_data.xlsx', sheet = 'info')
X2 <- rio::import('data/cell_mono_bi_data.xlsx', sheet = 'matrix') %>% 
  column_to_rownames('rownames')
info2 <- rio::import('data/cell_mono_bi_data.xlsx', sheet = 'info')


# histone
dfprot_his <- dfprot %>%
  filter(sapply(Protein.Group, function(x){
    any(unlist(str_split(x, ';')) %in% mouse_histone$Entry)
  }))
X_his <- df_faxp[dfprot_his$Protein.Group, ]
X_his_imp <- X_his %>%
  impute_min(., minvalue = apply(., 1, min, na.rm = T),
             downfactor = 0.9,
             SD = apply(., 1, sd, na.rm = T), width = 0.3)

df_his_detect <- X_his %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('ID') %>%
  pivot_longer(cols = -ID, names_to = 'Protein.Group', values_to = 'Log2Intensity') %>%
  mutate(Is.Detected = !is.na(Log2Intensity)) %>% 
  inner_join(dfprot, .) %>%
  inner_join(info, .) %>% 
  distinct(ID, Protein.Group, Is.Detected)
df_his <- X_his_imp %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('ID') %>%
  pivot_longer(cols = -ID, names_to = 'Protein.Group', values_to = 'Log2Intensity') %>%
  inner_join(dfprot, .) %>%
  inner_join(info, .) %>% 
  inner_join(., df_his_detect) %>% 
  mutate(Type = as.character(Type))

df_his[df_his$Type == '1nucleus', 'Type'] <- 'Nucleus'
df_his[df_his$Type == 'cell_1nucleus', 'Type'] <- 'Cell (Mononuclear)'
df_his[df_his$Type == 'cell_2nucleus', 'Type'] <- 'Cell (Binuclear)'
df_his$Type %<>% factor(levels = c('Nucleus', 'Cell (Mononuclear)', 'Cell (Binuclear)'))

comparisons <- combn(levels(df_his$Type), 2, simplify = F)


# visualization
df_his$Gene <- sapply(df_his$Genes, function(x){
  str_split(x, ';')[[1]][1]
})
set.seed(10)
pS8B <- ggplot(df_his) +
  aes(x = Gene, y = Log2Intensity, color = Type) +
  geom_boxplot(fill = '#FFFFFF', width = 0.8, outlier.shape = NA, position = position_dodge())+
  geom_point(aes(shape = Is.Detected, group = Type), alpha = 1, size = 3, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8)
  ) +
  scale_color_manual(values = type_colors) +
  scale_shape_manual(values = c('TRUE' = 20, 'FALSE' = 4)) +
  theme_classic() +
  theme(text = element_text(size = 15))
ggsave('Figure_S8B.pdf', pS8B, width = 12, height = 5)

tblS8B <- df_his


# Part II output -------
list(`Figure4F_#proteins` = tbl4F,
     Figure4G_pearson = tbl4G,
     Figure4H_GSEA_GOCC = tbl4H1,
     Figure4H_GSEA_Reactome = tbl4H2) %>% 
  rio::export('F4_FGH_source_data.xlsx')


p4F <- ggpubr::ggarrange(p4F, labels = 'F')
p4H <- ggpubr::ggarrange(p4H1, p4H2, common.legend = T, nrow = 1, ncol = 2, labels = 'H', legend = 'right')
ggsave('F4F_proteins.pdf', p4F, width = 5, height = 5)
ggsave('F4H_GSEA.pdf', p4H, width = 15, height = 8)



pS7 <- ggpubr::ggarrange(pS7A, pS7B, labels = c('A', 'B'))
ggsave('Figure_S7AB.pdf', pS7, width = 10, height = 5)
pS8B <- ggpubr::ggarrange(pS8B, labels = c('B'))
ggsave('Figure_S8B.pdf', pS8B, width = 5, height = 5)

list(FigureS7A = tblS7A,
     FigureS7B = tblS7B,
     FigureS7C = tblS7C,
     FigureS7D = tblS7D,
     FigureS7E = tblS7E,
     FigureS8A = tblS8A,
     FigureS8B_histone = tblS8B,
     FigureS8C_upset = tblS8C) %>% 
  rio::export('Figure_S7_S8_source_data.xlsx')


save.image('Figure4FGH_S7S8.RData')

