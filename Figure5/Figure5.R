source('../universal_functions/FAXP_init.R')

# 1.Figure D, F, G --------------------------------------------------------
## 1.1 load data -------
#protein matrix
dfpg <- rio::import('data/whole_library/Whole_lib_report.pg_matrix.tsv')
dfpr <- rio::import('data/whole_library/Whole_lib_report.pr_matrix.tsv')
colnames(dfpg) %<>% str_remove('^.+\\\\')
colnames(dfpr) %<>% str_remove('^.+\\\\')

#protein info
protinfo <- dfpg %>% select(Protein.Group:First.Protein.Description)

#protein matrix (log2 transform)
mat <- dfpg %>%
  select(-(Protein.Ids:First.Protein.Description)) %>%
  column_to_rownames('Protein.Group') %>% log2()

#sample info
info <- rio::import('data/20230219_Batch_design_CRC_demo.xlsx', sheet = 'Sheet2') %>% select(-ncol(.))
info <- data.frame(SampleName = colnames(mat),
                   Batch = str_extract(colnames(mat), 'b\\d+_(\\d+|pool)')) %>% 
  full_join(info)
info$BatchHead <- str_remove(info$Batch, '_.+$')
info$Slide %<>% factor()
info$Region <- ifelse(str_detect(info$Region, 'PC|CC'), 'C', info$Region)
info$Region %<>% factor(levels = region_order, ordered = T)

#rename files in mat and regenerate dfpg
colnames(mat) %<>% str_extract('b\\d+_(\\d+|pool)')
df <- mat %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>% 
  inner_join(info, .)


## 1.2 pool correlation ----------
pool_data <- list(
  pg = dfpg %>%
    select(-(Protein.Ids:First.Protein.Description)) %>%
    column_to_rownames('Protein.Group') %>% log2() %>% 
    select(matches('pool')) %>% 
    na_ratio_cutoff(cutoff = 1, MARGIN = 1) %>% 
    set_colnames(., str_extract(colnames(.), 'b\\d_pool')),
  pr = dfpr %>% 
    select(-(Protein.Group:Precursor.Charge)) %>%
    column_to_rownames('Precursor.Id') %>% log2() %>% 
    select(matches('pool')) %>% 
    na_ratio_cutoff(cutoff = 1, MARGIN = 1) %>% 
    set_colnames(., str_extract(colnames(.), 'b\\d_pool'))
)

cors_pg <- cor(pool_data$pg, use = 'pairwise.complete.obs', method = 'pearson')
cors_pr <- cor(pool_data$pr, use = 'pairwise.complete.obs', method = 'pearson')
pdf('Figure_S9C.pdf', width = 8, height = 8)
my_corrplot(cors_pg)
my_corrplot(cors_pr)
graphics.off()

tblS9C1 <- cors_pg %>% as.data.frame() %>% rownames_to_column('ID')
tblS9C2 <- cors_pr %>% as.data.frame() %>% rownames_to_column('ID')



## 1.3 protein group numbers (D) ----
### boxplot ----
dfbar <- apply(mat, 2, function(y) sum(!is.na(y))) %>% as.data.frame() %>% 
  setNames('# proteins') %>% 
  rownames_to_column('Batch') %>% 
  inner_join(info, .) %>%
  filter(!is.na(Patient)) # remove pooled samples

#remove outliers in each region by their # proteins
dfbar %<>% group_by(Region) %>% 
  summarise(lower = get_outliers(`# proteins`)[1],
            higher = get_outliers(`# proteins`)[2]) %>% 
  left_join(dfbar, .) %>% 
  mutate(Included = between(`# proteins`, lower, higher))

df %<>% filter(Batch %in% dfbar$Batch[dfbar$Included])
mat <- mat[, df$Batch]

dfbar.stat <- create_stat_value(dfbar, 'Region', '# proteins', padj.method = 'BH', nSignif = 2, test.fun = t.test, var.equal = T)

set.seed(10)
p5D <- ggplot(dfbar, aes(x = Region, y = `# proteins`, color = Region)) +
  stat_boxplot(geom = 'errorbar', width = 0.4)+
  geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
  geom_jitter(alpha = 0.7, size = 3, width = 0.35)+
  labs(x = '', y = '# proteins', subtitle = "Welch's t-test (BH-adjusted)")+
  scale_y_continuous(limits = c(0, 8200)) +
  scale_color_manual(values = region_colors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_pvalue_manual(
    data = dfbar.stat, label = 'p.adj'
  )
tbl5D <- dfbar %>% cbind(NA) %>% cbind(NA) %>% 
  cbind(rbind(dfbar.stat, matrix(rep(NA, (nrow(dfbar) - nrow(dfbar.stat)) * ncol(dfbar.stat)), ncol = ncol(dfbar.stat)) %>% set_colnames(colnames(dfbar.stat))))
names(tbl5D)[names(tbl5D) == 'NA'] <- ''


### heatmap ----
dfheat <- dfbar %>%
  filter(Included) %>% 
  select(Batch, Region, `# proteins`, Slide, Patient, BatchHead) %>%
  arrange(Region) %>%
  left_join(mat %>% t() %>% as.data.frame() %>% rownames_to_column('Batch')) %>% 
  set_rownames(.$Batch)
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

anno_colors <- list(
  Region = region_colors,
  BatchHead = brewer.pal(length(unique(dfheat$BatchHead)), 'Set2') %>% setNames(sort(unique(dfheat$BatchHead))),
  Patient = patient_colors,
  '# proteins' = rev(brewer.pal(11, 'RdBu')),
  Slide = brewer.pal(nlevels(dfheat$Slide), 'Oranges') %>% setNames(levels(dfheat$Slide))
)

clust_col <- hclust(dist(t(mat_heat), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()


bk <- unique(c(seq(-2, 2, length = 50)))
p5E <- pheatmap(mat_heat,
         annotation_col = dfheat %>% select(`# proteins`:BatchHead, Region),
         annotation_colors = anno_colors,
         cluster_rows = T, cluster_cols = clust_col,
         clustering_distance_rows = 'euclidean',
         clustering_distance_cols = 'euclidean',
         clustering_method = 'ward.D2',
         show_rownames = F, show_colnames = T,
         fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
         color = heat_colors, breaks = bk, border_color = F,
         filename = 'F5E_whole_heatmap.pdf',
         width = 10, height = 10
)
tbl5E <- mat_heat[p5E$tree_row$labels[p5E$tree_row$order], p5E$tree_col$labels[p5E$tree_col$order]] %>% as.data.frame() %>% rownames_to_column('Protein')


## 1.4 DEA (F) ----
### missing ratio control -----
# single protein level control: <20% in at least 3 Regions
df %>% count(Region)

df_miss <- df %>% group_by(Region) %>% select(-(SampleName:BatchHead)) %>%
  summarise_all(function(y) sum(is.na(y)) / length(y)) %>% 
  mutate(Region = as.character(Region))
df_miss %<>% rbind(c(Region = 'ALL', apply(mat, 1, function(x) sum(is.na(x)) / length(x))))

tmp <- df_miss %>%
  filter(Region != 'ALL') %>% 
  column_to_rownames('Region') %>%
  apply(2, function(y) sum(y < 0.2) >= 3 )
selected_pro <- names(tmp[tmp])
df_miss %<>% select(Region, all_of(selected_pro))
mat <- mat[colnames(df_miss)[-1], ]



# whole matrix level control
mr <- apply(mat, 1, function(x) sum(is.na(x)) / length(x))
plot(density(mr))

x <- seq(0.01, 0.99, 0.01)
y <- sapply(x, function(mr_cutoff){
  X <- mat[mr < mr_cutoff, ]
  return(sum(is.na(X)) / nrow(X) / ncol(X))
})
data.frame(x, y)
plot(x * 100, y * 100, xlab = '% NA ratio cutoff on protein level', ylab = '% whole matrix NA ratio')
abline(h = 6.5, v = 50)
# whole matrix NA ratio is 18.9% while the cutoff on protein level is set to 70%

mat50 <- mat[mr < 0.5, ]
dim(mat) # 3464  123
dim(mat50) # 3464  123


# NA ratio was considered by region respectively, and NA will be omitted
dfmiss <- df %>%
  filter(!is.na(Region)) %>% 
  group_by(Region) %>% 
  select(-(SampleName:BatchHead)) %>% 
  summarise_all(function(y) sum(is.na(y)) / length(y))

df_dea <- mat %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>%
  inner_join(info, .) %>%
  filter(!is.na(Region))


### limma ------------
library(limma)
sample_info <- df_dea %>% select(Batch, Patient, Slide, Region)
protein_expression_matrix <- df_dea %>%
  column_to_rownames('Batch') %>%
  select(-(SampleName:BatchHead)) %>% t()

comparison <- rev(c('RegionC-RegionH', 'RegionH-RegionL', 'RegionL-RegionN', 'RegionH-RegionN', 'RegionC-RegionL', 'RegionC-RegionN'))
design <- model.matrix(~0 + Region, data = sample_info) # Region as fixed effect
contrast <- makeContrasts(contrasts = comparison, levels = design)

dupcor <- duplicateCorrelation(protein_expression_matrix, design = design, block = sample_info$Patient)
dupcor$consensus.correlation # 0.1518476
plot(density(dupcor$atanh.correlations, bw = 'nrd0'), main = 'atanh.correlations', xlab = '\n(range [-0.0028, 1.2364] 15% trimmed mean = 0.1518)\n(N = 3464 Bandwidth = 0.03329)')
abline(v = dupcor$consensus.correlation)
abline(v = min(dupcor$atanh.correlations))
abline(v = max(dupcor$atanh.correlations))

# Estimate the treatment effect using both complete and incomplete blocks
fit1 <- lmFit(protein_expression_matrix, design = design, block = sample_info$Patient, correlation = dupcor$consensus)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
res_all <- topTable(fit3, p.value = 0.05, adjust = 'fdr', sort.by = 'B', number = Inf) %>% 
  rownames_to_column('Protein')

tbl6E <- res_all %>% select(Protein, matches('\\.RegionN$')) %>%
  pivot_longer(cols = -Protein, names_to = 'Comparison', values_to = 'Log2FC') %>% 
  filter(abs(Log2FC) > log2(2)) %>% 
  inner_join(res_all %>% distinct(Protein, adj.P.Val)) %>% 
  inner_join(protinfo %>% rename(Protein = Protein.Ids) %>% distinct(Protein, Genes), .) %>% 
  arrange(Comparison, desc(Log2FC))

res_ls <- lapply(setNames(colnames(contrast), colnames(contrast)), function(coef){
  topTable(fit3, coef = coef, p.value = 1, adjust = 'fdr', sort.by = 'B', number = Inf) %>% rownames_to_column('Protein')
})


#### missing value filter -----
df_quant <- lapply(comparison, function(comparex){
  comp <- comparex %>% str_remove_all('Region') %>% str_split('-') %>% .[[1]]
  X <- df_dea %>% filter(Region %in% comp) %>%
    select(-(SampleName:Slide), -(Rep:BatchHead)) %>% 
    pivot_longer(cols = -Region, names_to = 'Protein', values_to = 'Log2Intensity')
  ret <- plyr::ddply(X, 'Protein', function(dfsub){
    x1 <- dfsub$Log2Intensity[dfsub$Region == comp[1]]
    x2 <- dfsub$Log2Intensity[dfsub$Region == comp[2]]
    # Log2FC <- log2(mean(2^x1, na.rm = T) / mean(2^x2, na.rm = T))
    # P <- tryCatch(t.test(x1, x2)$p.value, error = function(e) NA)
    ret <- data.frame(Comparison = comparex#,
                      # Log2FC = Log2FC,
                      #P = P
    )
    # ret$adj.P <- p.adjust(ret$P, 'BH')
    return(ret)
  })
}) %>% plyr::ldply()

DEA <- protinfo %>%
  select(Protein.Ids, Genes) %>%
  rename(Protein = Protein.Ids) %>%
  inner_join(df_quant) %>%
  inner_join(dfmiss %>% column_to_rownames('Region') %>% t() %>% as.data.frame() %>% setNames(str_c('MissingRatio_', names(.))) %>% rownames_to_column('Protein')) %>% 
  inner_join(plyr::ldply(res_ls, .id = 'Comparison')) %>% 
  inner_join(
    res_all %>%
      select(Protein, AveExpr:adj.P.Val) %>%
      rename_at(vars(-1), function(x) str_c(x, '.Global'))
  )
DEA %<>% 
  rename(Log2FC = logFC,
         P = P.Value,
         adj.P = adj.P.Val) %>% 
  mutate(
    Type = ifelse(Log2FC > 0, 'Up', 'Down'),
    is.significant = abs(Log2FC) > log2(2) & adj.P < 0.05 & adj.P.Val.Global < 0.05,
    Label = str_c(Genes, ' (', Protein, ')'),
    Comparison = factor(Comparison, levels = comparison, ordered = T)
  ) %>% 
  arrange(Comparison) %>% 
  mutate(Comparison = factor(
    str_remove_all(Comparison, 'Region'),
    levels = str_remove_all(levels(Comparison), 'Region')
  ))

DEP <- DEA %>% filter(is.significant)



## 1.4.x DEA visualization -----
#p value profile
ggplot(DEA) +
  aes(x = Comparison, y = adj.P, fill = Comparison) +
  geom_violin(alpha = 0.7, width = 0.5) +
  geom_hline(yintercept = 0.05, color = 'red') +
  scale_y_log10() +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic()

### DEP number ------
tmp <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)[c(1,2,6)]) %>%
  mutate(Comparison = factor(Comparison, levels = rev(levels(DEA$Comparison)))) %>%
  select(Label, Comparison, Log2FC, adj.P, Type, is.significant) %>%
  count(Comparison, Type)
tbl6C_line1 <- dfline1 <- tmp %>% group_by(Comparison) %>% summarise(nn = sum(n)) %>% 
  inner_join(tmp, .) %>% 
  mutate(x = str_c(Comparison, '\n(', nn, ')')) %>% 
  arrange(Comparison) %>% 
  mutate(x = factor(x, unique(x)))
p6C1 <- ggplot(dfline1, aes(x = x, y = n, group = Type, color = Type))+
  geom_line(linewidth = 1) +
  geom_point(size = 4, alpha = 1) +
  labs(x = '', y = '', subtitle = '# differentially expressed proteins')+
  scale_color_manual(values = c(Up = 'red4', Down = 'blue4')) +
  scale_y_continuous(limits = c(0, 150)) +
  theme_classic() +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 0, hjust = 0.5))

tmp <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)[c(1,3,4)]) %>%
  mutate(Comparison = factor(Comparison, levels = rev(levels(DEA$Comparison)))) %>%
  select(Label, Comparison, Log2FC, adj.P, Type, is.significant) %>%
  count(Comparison, Type)
tbl6C_line2 <- dfline2 <- tmp %>% group_by(Comparison) %>% summarise(nn = sum(n)) %>% 
  inner_join(tmp, .) %>% 
  mutate(x = str_c(Comparison, '\n(', nn, ')')) %>% 
  arrange(Comparison) %>% 
  mutate(x = factor(x, rev(unique(x))))
p6C2 <- ggplot(dfline2, aes(x = x, y = n, group = Type, color = Type))+
  geom_line(linewidth = 1) +
  geom_point(size = 4, alpha = 1) +
  labs(x = '', y = '', subtitle = '# differentially expressed proteins')+
  scale_color_manual(values = c(Up = 'red4', Down = 'blue4')) +
  theme_classic() +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 0, hjust = 0.5))

# tbl6C_line <- dfline <- DEP %>%
#   filter(Comparison %in% levels(DEA$Comparison)) %>%
#   select(Label, Comparison, Log2FC, adj.P, Type, is.significant) %>%
#   count(Comparison, Type)



### heatmap -------
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

heat1 <- pheatmap(mat_heat[unique(DEP$Protein), ],
                  annotation_col = dfheat %>% select(Region),
                  annotation_colors = anno_colors,
                  cluster_rows = T, cluster_cols = F,
                  clustering_distance_rows = 'euclidean',
                  clustering_distance_cols = 'euclidean',
                  clustering_method = 'ward.D2',
                  show_rownames = F, show_colnames = F
)

mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score


clust_col <- hclust(dist(t(mat_heat[heat1$tree_row$labels[heat1$tree_row$order], ]), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()


p6A <- pheatmap(
  mat_heat[heat1$tree_row$labels[heat1$tree_row$order], ],
  scale = 'none',
  annotation_col = dfheat %>% select(Region, Patient),
  annotation_colors = anno_colors,
  cluster_rows = T, cluster_cols = clust_col,
  clustering_method = 'ward.D2',
  clustering_distance_cols = 'euclidean', clustering_distance_rows = 'euclidean',
  show_rownames = F, show_colnames = F,
  color = heat_colors, breaks = bk,
  fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
  border_color = F, na_col = '#AAAAAA',
  filename = 'F6A_335DEPs_heatmap_hclust.pdf',
  width = 10, height = 5)
tbl6A <- mat_heat[p6A$tree_row$labels[p6A$tree_row$order], p6A$tree_col$labels[p6A$tree_col$order]] %>% as.data.frame() %>% rownames_to_column('Protein')


### UMAP ------
library(umap)
X_imp <- mat_heat[heat1$tree_row$labels[heat1$tree_row$order], ]
X_imp[is.na(X_imp)] <- 0

umap <- umap(t(X_imp), n_neighbors = 56, min_dist = 0.5) # default n_neighbors = 10, min_dist = 0.1
df_umap <- info %>% inner_join(umap$layout %>% as.data.frame() %>% rownames_to_column('Batch'))
df_umap_center <- df_umap %>% group_by(Region) %>% summarise_at(vars(center_x = V1, center_y = V2), mean) # calculate center of clusters
df_umap %<>% left_join(df_umap_center)

p6B <- ggplot(df_umap, aes(x = V1, y = V2, color = Region)) +
  geom_point(size = 4, alpha = 0.9)+
  labs(x = 'UMAP1', y = 'UMAP2') +
  scale_color_manual(values = region_colors) +
  theme_classic() +
  theme(text = element_text(size = 10), legend.position = 'none')

ggplot(df_umap, aes(x = V1, y = V2, color = Patient)) +
  geom_point(size = 4, alpha = 0.9)+
  labs(x = 'UMAP1', y = 'UMAP2') +
  scale_color_manual(values = patient_colors) +
  theme_classic() +
  theme(text = element_text(size = 10), legend.position = 'right')

tbl6B <- df_umap


### volcano -----
set.seed(10)
p6C_volcano1 <- DEA %>% 
  filter(Comparison %in% levels(DEA$Comparison)[c(1,2,6)]) %>% 
  mutate(Comparison = factor(Comparison, rev(levels(DEA$Comparison)))) %>% 
  ggplot(aes(x = Comparison, y = Log2FC, color = interaction(Type, is.significant)))+
  geom_jitter(alpha = 0.6, size = 3, width = 0.4) +
  labs(x = '', y = 'Log2 Fold-Change', subtitle = '# differentially expressed proteins')+
  scale_color_manual(name = '',
                     values = c(Up.FALSE = '#AAAAAA', Down.FALSE = '#AAAAAA',
                                Up.TRUE = 'red4', Down.TRUE = 'blue4'),
                     labels = c('', '|FC|<2 or adj.P>0.05', 'Down', 'Up')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'top') #+

set.seed(10)
p6C_volcano2 <- DEA %>% 
  filter(Comparison %in% levels(DEA$Comparison)[c(1,3,4)]) %>% 
  ggplot(aes(x = Comparison, y = Log2FC, color = interaction(Type, is.significant)))+
  geom_jitter(alpha = 0.6, size = 3, width = 0.4) +
  labs(x = '', y = 'Log2 Fold-Change', subtitle = '# differentially expressed proteins')+
  scale_color_manual(name = '',
                     values = c(Up.FALSE = '#AAAAAA', Down.FALSE = '#AAAAAA',
                                Up.TRUE = 'red4', Down.TRUE = 'blue4'),
                     labels = c('', '|FC|<2 or adj.P>0.05', 'Down', 'Up')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'top') #+

tbl6C <- DEA


# 2.IPA analysis ----------------------------------------------------------
cal_med_pnts <- function(x){
  # Calculate mid points of x, after add 0 at the begin
  # x is an an atomic vector;
  x <- c(0, x)
  sapply(2:length(x), function(i){
    (x[i-1] + x[i]) / 2
  })
}


## Biomarker
tbls <- list.files('data/IPA_analysis/', pattern = '\\.xls$', full.names = T) %>%
  setNames(., str_extract(., '/([A-Za-z]+vs[A-Za-z]+)', group = 1))

ipa <- plyr::ldply(tbls, function(tbl){
  df <- rio::import(tbl)
  colnames(df) <- df[1, ]
  df <- df[-1, ]
  df$`Expr Fold Change` %<>% as.numeric()
  return(df)
})
df_ipa <- ipa %>% select(.id, Symbol, `Biomarker Application(s)`) %>% 
  rename(Biomarker = `Biomarker Application(s)`) %>% 
  separate_rows('Biomarker', sep = ',') %>% 
  drop_na()
df_ipa$Biomarker <- factor(str_to_sentence(df_ipa$Biomarker), levels = c('Efficacy', 'Response to therapy', 'Safety', 'Diagnosis', 'Disease progression', 'Prognosis', 'Unspecified application'))
df_ipa$.id <- factor(df_ipa$.id, levels = c('CCvsN', 'HvsN', 'LvsN'))
df_ipa_stat <- df_ipa %>%
  count(Biomarker, .id)

p6D <- ggplot(df_ipa_stat, aes(x = Biomarker, y = n, fill = .id)) +
  geom_col() +
  geom_text(aes(label = n), size = 3,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = region_colors[c('L', 'H', 'C')] %>% setNames(c('LvsN', 'HvsN', 'CCvsN'))) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Biomarker Application', y = 'Count') +
  theme_classic() +
  theme(text = element_text(size = 10), line = element_line(linewidth = 0.5), legend.position = 'right')


tbl6D <- df_ipa %>% inner_join(df_ipa_stat)





# output ------
library(patchwork)
list(Figure5D_protein = tbl5D,
     Figure5E_heatmap = tbl5E,
     Figure6A_heatmapDEP = tbl6A,
     Figure6B_UMAP = tbl6B,
     `Figure6C_DEA(&SF10A)` = tbl6C,
     Figure6D_IPA = tbl6D,
     Figure6E_marker = tbl6E) %>% 
  rio::export('Figure_5DE_6A-E_source_data.xlsx')


p5D <- ggpubr::ggarrange(p5D, nrow = 1, labels = 'D')
ggsave('F5D_protein_number.pdf', p5D, width = 5, height = 5)  

p6BD <- ggpubr::ggarrange(p6B, p6D, nrow = 1, labels = c('B', 'D'))
ggsave('F6BD.pdf', p6BD, width = 10, height = 5)

p6C_line <- ggpubr::ggarrange(p6C1, p6C2, nrow = 1, ncol = 2, common.legend = T)
p6C_volcano <- ggpubr::ggarrange(p6C_volcano1, p6C_volcano2, nrow = 1, ncol = 2, common.legend = T)
p6C <- wrap_plots(p6C_volcano, p6C_line) + plot_layout(heights = c(3.5,1))
ggsave('Figure_6C_S10A.pdf', p6C, width = 12, height = 7.7)




list(FigureS9C_protein = tblS9C1,
     FigureS9C_peptide = tblS9C2) %>% 
  rio::export('Figure_S9C_source_data.xlsx')



save.image('Figure_5DE_6A-E_S9C_S10A.RData')



