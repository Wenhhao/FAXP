source('../universal_functions/FAXP_init.R')


# 1.B-J -------------------------------------------------------------------
## 1.1 load data ----
#peptide matrix
dfpep <- rio::import('data/InTip_vs_InGel_DDA/combined_peptide.tsv')
dfpep %<>% filter(!str_detect(Protein, '^CON_'))

matpep <- dfpep %>%
  select(`Peptide Sequence`, matches('^In\\w+?_1_\\d Intensity$')) %>%
  column_to_rownames('Peptide Sequence')
colnames(matpep) %<>% str_remove(' Intensity$')

matpep_count <- dfpep %>%
  select(`Peptide Sequence`, matches('^In\\w+?_1_\\d.+Spectral Count$')) %>%
  column_to_rownames('Peptide Sequence')
colnames(matpep_count) %<>% str_remove(' Spectral Count$')

#sample info
info <- rio::import('data/F3_info.xlsx', sheet = 1)
info <- data.frame(SampleName = colnames(matpep),
                   'Label in file' = str_remove(colnames(matpep), '\\d+$'),
                   check.names = F) %>% left_join(info)

#protein matrix
dfpro <- rio::import('data/InTip_vs_InGel_DDA/combined_protein.tsv')
dfpro %<>% filter(!str_detect(Protein, '^CON_'))
matpro <- dfpro %>%
  mutate(tmp = `Indistinguishable Proteins` %>%
           sapply(function(x){
             str_split(x, ', ')[[1]] %>%
               str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
               str_c(collapse = ', ')
           }),
         `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
  select(`Protein ID`, matches('^In\\w+?_1_\\d Intensity$')) %>%
  column_to_rownames('Protein ID')
colnames(matpro) %<>% str_remove(' Intensity$')

matpro_count <- dfpro %>%
  mutate(tmp = `Indistinguishable Proteins` %>%
           sapply(function(x){
             str_split(x, ', ')[[1]] %>%
               str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
               str_c(collapse = ', ')
           }),
         `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
  select(`Protein ID`, matches('^In\\w+?_1_\\d Total Spectral Count$')) %>%
  column_to_rownames('Protein ID')
colnames(matpro_count) %<>% str_remove(' Total Spectral Count$')

#regenerate dfpep and dfpro; within protein groups
dfpep <- matpep %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro <- matpro %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpep_count <- matpep_count %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro_count <- matpro_count %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
pep_list <- plyr::dlply(dfpep_count, '`Label in figure`', function(dfsub){
  dfsub %>% select(-(SampleName:Note)) %>% 
    .[, apply(., 2, function(y) !all(y == 0))] %>% 
    colnames()
})
pro_list <- plyr::dlply(dfpro_count, '`Label in figure`', function(dfsub){
  dfsub %>% select(-(SampleName:Note)) %>% 
    .[, apply(., 2, function(y) !all(y == 0))] %>% 
    colnames() %>% str_split(', ') %>% unlist() %>% unique()
})


## 1.2 identity (B) ---------
### 1.2.1 peptide number (B) ---------
dflong <- dfpep_count %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Peptide', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA))
tbl1a <- dfbar <- dflong %>% 
  group_by(SampleName) %>% 
  summarise(`# Peptides` = sum(!is.na(Count))) %>% 
  left_join(info)

set.seed(0)
p1a <- ggplot(dfbar, aes(x = `Label in figure`, y = `# Peptides`, color = `Label in figure`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '# peptides')+
  scale_color_manual(values = gel_tip_colors) +
  scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)

### 1.2.2 protein number (B) ---------
dflong <- dfpro_count %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Protein', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA)) %>% 
  separate_rows('Protein')
tbl1b <- dfbar <- dflong %>% 
  group_by(SampleName) %>% 
  summarise(`# Proteins` = sum(!is.na(Count))) %>% 
  left_join(info)

set.seed(0)
p1b <- ggplot(dfbar, aes(x = `Label in figure`, y = `# Proteins`, color = `Label in figure`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '# proteins')+
  scale_color_manual(values = gel_tip_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)

tbl1 <- tbl1a %>% inner_join(tbl1b) %>% select(1, 3, 4, 2, 6)
p1 <- ggpubr::ggarrange(p1a, p1b)



## 1.3 Missed cleavage (%) (C) --------
pepList <- dfpep_count %>% 
  select(-Note) %>% 
  plyr::dlply('SampleName', function(dfsub){
    dfsub[, which(dfsub != 0)] %>% select(-(SampleName:'Label in figure')) %>% colnames()
  })
df_missCleav <- plyr::ldply(pepList, calc_missCleav, .id = 'SampleName')
tbl2 <- dfbar <- df_missCleav %>%
  count(SampleName, MissedCleavage) %>%
  with_groups(SampleName, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  inner_join(info)

set.seed(0)
p2 <- ggplot(dfbar, aes(x = `Label in figure`, y = `ratio (%)`, color = `Label in figure`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage)), scales = 'free') +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  scale_color_manual(values = gel_tip_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF')) +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)


## 1.4 Physicochemical property (D) --------
dfscore1 <- plyr::ldply(pep_list, calc_pep1, .id = 'Label in figure')
dfscore2 <- plyr::ldply(pep_list, calc_pep2, .id = 'Label in figure')
tbl3 <- dfscore <- dfscore1 %>% inner_join(dfscore2)

p3_list <- list()
for(i in 3:ncol(dfscore)){
  dftmp <- dfscore %>% select(1, all_of(i)) %>% setNames(c('Label in figure', 'Index'))
  p3_list[[i-2]] <- ggplot(dftmp, aes(x = Index, color = `Label in figure`)) +
    geom_density(alpha = 0.9, linewidth = 1) +
    labs(x = '', y = '', subtitle = str_c('Density of ', colnames(dfscore)[i]))+
    scale_color_manual(values = gel_tip_colors) +
    theme_classic() +
    theme(text = element_text(size = 15), legend.position = 'none')
}
names(p3_list) <- colnames(dfscore)[-(1:2)]
p3_all <- ggpubr::ggarrange(plotlist = p3_list, nrow = 5, ncol = 5)
ggsave('F3D_peptide_physicochemical_property.pdf', p3_all, width = 14, height = 14)
p3 <- p3_list$pI_EMBOSS


## 1.5 Peptide length curve (E) --------
df <- rio::import('data/InTip_vs_InGel_DDA/combined_peptide.tsv')
df %<>% filter(!str_detect(`Protein`, '^CON_')) %>% select(`Peptide Sequence`, `Peptide Length`) %>% rename(Peptide = `Peptide Sequence`)

dflong <- dfpep_count %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Peptide', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA))
tbl4 <- dfbar <- dflong %>% 
  left_join(df) %>% 
  filter(!is.na(Count)) %>% 
  select(-SampleName, -Note, -Count) %>% 
  distinct() %>% 
  count(`Label in figure`, `Peptide Length`)

p4 <- ggplot(dfbar, aes(x = `Peptide Length`, y = n,
                        color = `Label in figure`)) +
  geom_point(size = 2) +
  labs(x = '', y = '', subtitle = 'Peptide length')+
  scale_y_continuous(labels = scales::scientific) +
  scale_color_manual(values = gel_tip_colors) +
  scale_fill_manual(values = gel_tip_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')


## 1.6 Venn (F) --------
# VennDiagram
venn1 <- VennDiagram::venn.diagram(x = pep_list,
                                resolution = 300,
                                alpha=rep(0.95, length(pep_list)),
                                # fill=allFills[c(1, 4, 5)],
                                fill = 'white', col = gel_tip_colors,
                                main = stringr::str_glue("Peptides ({nrow(matpep)} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
venn2 <- VennDiagram::venn.diagram(x = pro_list,
                                resolution = 300,
                                alpha=rep(0.95, length(pro_list)),
                                fill = 'white', col = gel_tip_colors,
                                main=stringr::str_glue("Proteins ({nrow(matpro)} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
pdf('F3F_VennDiagram.pdf', width = 10, height = 10)
grid::grid.newpage(); grid::grid.draw(venn1)
grid::grid.newpage(); grid::grid.draw(venn2)
graphics.off()


vennlist <- list(InGel.Peptide = pep_list$InGel,
                 InTip.Peptide = pep_list$InTip,
                 InGel.Protein = pro_list$InGel,
                 InTip.Protein = pro_list$InTip)
venntbl <- matrix(nrow = max(sapply(vennlist, length)), ncol = 4) %>% 
  as.data.frame() %>% setNames(names(vennlist))
for(j in 1:ncol(venntbl)){
  venntbl[1:length(vennlist[[j]]), j] <- vennlist[[j]]
}
tbl5 <- venntbl



## 1.7 CV (G) --------
### protein CV -----
dfcv_pro <- dfpro %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Protein', values_to = 'Intensity') %>% 
  mutate(Intensity = ifelse(Intensity != 0, Intensity, NA)) %>% 
  filter(!is.na(Intensity)) %>% 
  group_by(`Label in figure`, Protein) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean)

plot_procv <- ggplot(dfcv_pro, aes(x = `Label in figure`, y = CV, color = `Label in figure`, fill = `Label in figure`))+
  geom_violin(color = NA, width = 1)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Proteins')+
  scale_color_manual(values = gel_tip_colors) +
  scale_fill_manual(values = gel_tip_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)
x1 <- dfcv_pro %>% filter(`Label in figure` == 'InGel') %>% pull(CV)
x2 <- dfcv_pro %>% filter(`Label in figure` == 'InTip') %>% pull(CV)
lbl <- t.test(x1, x2)$p.value # 6.23e-43

### peptide CV -----
dfcv_pep <- dfpep %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Peptide', values_to = 'Intensity') %>% 
  mutate(Intensity = ifelse(Intensity != 0, Intensity, NA)) %>% 
  filter(!is.na(Intensity)) %>% 
  group_by(`Label in figure`, Peptide) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean)

plot_pepcv <- ggplot(dfcv_pep, aes(x = `Label in figure`, y = CV, color = `Label in figure`, fill = `Label in figure`))+
  geom_violin(color = NA, width = 1)+
  # stat_boxplot(geom = 'errorbar', color = '#000000', width = 0.05)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Peptides')+
  scale_color_manual(values = gel_tip_colors) +
  scale_fill_manual(values = gel_tip_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)
x1 <- dfcv_pep %>% filter(`Label in figure` == 'InGel') %>% pull(CV)
x2 <- dfcv_pep %>% filter(`Label in figure` == 'InTip') %>% pull(CV)
lbl <- t.test(x1, x2)$p.value # 0

p6 <- ggpubr::ggarrange(plot_pepcv, plot_procv)

tbl6a <- dfcv_pep %>%
  mutate(Molecular = 'Peptide', .before = 2) %>%
  rename(ID = Peptide)
tbl6b <- dfcv_pro %>%
  mutate(Molecular = 'Protein', .before = 2) %>%
  rename(ID = Protein)
tbl6 <- rbind(tbl6a, tbl6b)
rm(tbl6a, tbl6b)


## 1.8 Pearson's correlation (H) -----
matpro_ <- matpro
matpro_[matpro_ == 0] <- NA
cors_pro <- matpro_ %>% log2() %>%
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>%
  round(2)
min(cors_pro) # 0.84

matpep_ <- matpep
matpep_[matpep_ == 0] <- NA
cors_pep <- matpep_ %>% log2() %>%
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>%
  round(2)
min(cors_pep) # 0.56


pdf('F3H_pearson_correlation.pdf', width = 8, height = 8)
my_corrplot(cors_pro)
my_corrplot(cors_pep)
graphics.off()

list(Figure3H_Protein = cors_pro %>% as.data.frame() %>% rownames_to_column(),
     Figure3H_Peptide = cors_pep %>% as.data.frame() %>% rownames_to_column()) %>% rio::export('F3H_pearson_correlation.xlsx')


## 1.9 Log2Intensity rank, protein group level (I) --------
#Protein info are downloaded from UniProt.Org, and selected with "liver" detected in "Tissue specificity"
protinfo <- rio::import('data/uniprotkb_mouse_reviewed_17184_20231213.xlsx') %>% 
  filter(Entry %in% unique(unlist(str_split(colnames(dfpro)[-(1:4)], ', ')))) %>%
  rename(Protein = Entry) %>%
  select(-Reviewed, -Organism)

liverprotinfo <- protinfo %>%
  filter(str_detect(str_to_lower(`Tissue specificity`), 'liver'))

dfrank <- dfpro %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Protein', values_to = 'Intensity') %>% 
  mutate(Intensity = ifelse(Intensity != 0, Intensity, NA)) %>% 
  filter(!is.na(Intensity)) %>% 
  group_by(`Label in figure`, Protein) %>% 
  summarise(MeanIntensity = mean(Intensity, na.rm = T)) %>% 
  mutate(Log2MeanIntensity = log2(MeanIntensity)) %>% 
  arrange(desc(Log2MeanIntensity)) %>% 
  with_groups(`Label in figure`, mutate, Rank = 1:length(MeanIntensity))

tbl8 <- dfrank_label <- dfrank %>%
  filter(str_detect(Protein, str_c(liverprotinfo$Protein, collapse = '|')),
         !str_detect(Protein, ', ')) %>%  # do not label protein groups
  left_join(liverprotinfo)

p8 <- ggplot(dfrank, aes(x = Rank, y = Log2MeanIntensity, color = `Label in figure`, fill = `Label in figure`))+
  geom_line(color = '#000000', linewidth = 1) +
  geom_point(data = dfrank_label, size = 1, alpha = 0.5) +
  annotate(geom = 'text', x = 2000, y = 25,
           label = 'Liver specific proteins in Uniprot by in-gel (n=486)\nLiver specific proteins in Uniprot by in-tip (n=616)') +
  labs(x = 'Protein abundance rank', y = 'Protein level (log2)', subtitle = 'Protein group')+
  scale_color_manual(values = gel_tip_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')


## 1.10 Protein types&locations (J) --------
library(circlize)


IPA_ingel <- rio::import('data/InGel_all_molecules.xls')
IPA_intip <- rio::import('data/InTip_all_molecules.xls')

colnames(IPA_ingel) <- IPA_ingel[1, ]
IPA_ingel <- IPA_ingel %>% slice(-1) %>% rename(Type = `Type(s)`)
colnames(IPA_intip) <- IPA_intip[1, ]
IPA_intip <- IPA_intip %>% slice(-1) %>% rename(Type = `Type(s)`)

### ChordDiagram ----------
# for InGel
mat_ingel <- IPA_ingel %>%
  mutate(Location = factor(Location, levels = c('Nucleus', 'Cytoplasm', 'Plasma Membrane', 'Extracellular Space', 'Other'))) %>% 
  mutate(Type = ifelse(Type == 'other', 'other types', Type)) %>% 
  reshape2::dcast(Location ~ Type) %>%
  column_to_rownames('Location') %>% 
  t()
rownames(mat_ingel) %<>% str_to_sentence()
colnames(mat_ingel) %<>% str_to_sentence()

grid.col.location <- brewer.pal(8, 'Set2')[1:5]
names(grid.col.location) <- colnames(mat_ingel)

# grid.col.type <- rep('grey', nrow(mat_ingel))
grid.col.type <- my_color_palette(nrow(mat_ingel), brewer.pal(8, 'Spectral'), visible = F)
names(grid.col.type) <- rownames(mat_ingel)
grid.col <- c(grid.col.location, grid.col.type)

# border_mat_ingel <- matrix('black', nrow = 1, ncol = ncol(mat_ingel))
# rownames(border_mat_ingel)

pdf('F3J_chordDiagram_InGel.pdf', width = 5, height = 5)
chordDiagram(mat_ingel, grid.col = grid.col,
             # annotationTrack = "grid",
             link.lwd = 2, # width
             link.lty = 2, # style
             transparency = 0.5,
             big.gap = 10,
             small.gap = 2,
             # link.border = border_mat_ingel,
)
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  # circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.7, niceFacing = TRUE)
}, bg.border = NA)
title(str_glue('IPA of InGel'), cex = 0.8)
graphics.off()
circos.clear()


# for InTip
mat_intip <- IPA_intip %>%
  mutate(Location = factor(Location, levels = c('Nucleus', 'Cytoplasm', 'Plasma Membrane', 'Extracellular Space', 'Other'))) %>% 
  mutate(Type = ifelse(Type == 'other', 'other types', Type)) %>% 
  reshape2::dcast(Location ~ Type) %>%
  column_to_rownames('Location') %>% 
  t()
rownames(mat_intip) %<>% str_to_sentence()
colnames(mat_intip) %<>% str_to_sentence()

grid.col.location <- brewer.pal(8, 'Set2')[1:5]
names(grid.col.location) <- colnames(mat_intip)

# grid.col.type <- rep('grey', nrow(mat_intip))
grid.col.type <- my_color_palette(nrow(mat_intip), brewer.pal(8, 'Spectral'), visible = F)
names(grid.col.type) <- rownames(mat_intip)
grid.col <- c(grid.col.location, grid.col.type)

# border_mat_intip <- matrix('black', nrow = 1, ncol = ncol(mat_intip))
# rownames(border_mat_intip)

pdf('F3J_chordDiagram_InTip.pdf', width = 5, height = 5)
chordDiagram(mat_intip, grid.col = grid.col,
             # annotationTrack = "grid",
             link.lwd = 2, # width
             link.lty = 2, # style
             transparency = 0.5,
             big.gap = 10,
             small.gap = 2,
             # link.border = border_mat_intip,
)
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  # circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.7, niceFacing = TRUE)
}, bg.border = NA)
title(str_glue('IPA of InTip'), cex = 0.8)
graphics.off()
circos.clear()

plot_legend1 <- data.frame(Location = colnames(mat_ingel), y = 1) %>% 
  ggplot()+
  geom_col(aes(x = Location, y = y, fill = Location), width = 1)+
  scale_fill_manual(values = grid.col.location)+
  coord_flip()+
  theme_void()
plot_legend2 <- data.frame(Type = rownames(mat_ingel), y = 1) %>% 
  ggplot()+
  geom_col(aes(x = Type, y = y, fill = Type), width = 1)+
  scale_fill_manual(values = grid.col.type)+
  coord_flip()+
  theme_void()
plot_legend <- ggpubr::ggarrange(plot_legend1, plot_legend2, ncol = 1)
ggsave('F3J_chordDiagram_legends.pdf', plot_legend, width = 8, height = 8)



## output tables 2B-2J, and other figures -------
tblBtoI <- list(`Figure3B_#peptides_proteins` = tbl1,
                Figure3C_missed_cleavage = tbl2,
                Figure3D_physicochemistry = tbl3,
                Figure3E_peptide_length = tbl4,
                Figure3F_Venn = tbl5,
                Figure3G_CV_peptides = dfcv_pep,
                Figure3G_CV_proteins = dfcv_pro,
                Figure3H_Pearson_peptides = cors_pep %>% as.data.frame() %>% rownames_to_column(),
                Figure3H_Pearson_proteins = cors_pro %>% as.data.frame() %>% rownames_to_column(),
                Figure3I_abundance_rank = tbl8,
                Figure3J_InGel_chord_diagram = IPA_ingel,
                Figure3J_InTip_chord_diagram = IPA_intip)
rio::export(tblBtoI, 'F3_B-J_source_data.xlsx')

# figures
p <- ggpubr::ggarrange(p1, p2, p3, p4, p6, p8,
                       labels = LETTERS[c(2, 3, 4, 5, 7, 9)],
                       nrow = 2, ncol = 3)
ggsave('F3_B-J_others.pdf', p, width = 15, height = 10)


save.image('F3_B-J.RData')



