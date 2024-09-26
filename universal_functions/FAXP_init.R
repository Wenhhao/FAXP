pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(corrplot)
library(pheatmap)
library(ggsci)


# universal values --------------------------------------------------------
nucleus_type_order <- c('1nucleus', 'cell_1nucleus', 'cell_2nucleus', 'blank_control')
region_order <- c('N', 'L', 'H', 'C')

project_colors <- c(ProteomEx = '#00008B', FAXP = '#C32022')
gel_tip_colors <- c(InGel = '#00008B', InTip = '#C32022')
cor_colors <- viridis::viridis(101, alpha = 1, begin = 0.2, end = 1, direction = -1, option = 'A') %>% setNames(seq(0, 1, length.out = 101))
detect_range_colors <- c('In-tip DDA-MS' = '#DF8F44FF', 'In-gel DDA-MS' =  '#008280FF', 'In-tip DIA-MS' = '#B24745FF', 'In-gel DIA-MS' = '#00A1D5FF')
type_colors <- c("#358DB9FF", "#CF4E9CFF", "#8C57A2FF", "#2E2A2BFF") %>% 
  setNames(nucleus_type_order)
type_colors['Intra.Group.Nucleus'] <- type_colors['1nucleus']
type_colors['Intra.Group.Mono'] <- type_colors['cell_1nucleus']
type_colors['Intra.Group.Bi'] <- type_colors['cell_2nucleus']
type_colors['Inter.Groups'] <- type_colors['blank_control']
type_colors['Nucleus'] <- type_colors['1nucleus']
type_colors['Cell (Mononuclear)'] <- type_colors['cell_1nucleus']
type_colors['Cell (Binuclear)'] <- type_colors['cell_2nucleus']
heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu")))(50)
region_colors <- c(N = '#000080',
                   L = '#FFD700',
                   H = '#CC5500',
                   C = '#C71585'#,
                   #PC = '#E85ABD',
                   #CC = '#751E5B'
)
patient_colors <- c('#77B5FE', '#4682B4', '#004488') %>% setNames(str_c('P', 1:3))
pca_colors <- c('#D7DA86', '#330A5FFF', '#FCB519FF', '#ED6925FF', '#781C6DFF', '#BB3754FF', '#000004FF') %>% setNames(c('resid', 'Region', 'Patient:Region', 'Patient', 'Slide:Region', 'Patient:Slide', 'Slide'))


# universal functions -----------------------------------------------------
cv <- function(x, na.rm = T){
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}

get_outliers <- function(vec, coef = 1.5){
  # outliers based on Q1 and Q3
  stats <- quantile(vec, na.rm = T)
  iqr <- diff(stats[c(2, 4)])
  ret <- c(stats[2] - coef * iqr, stats[4] + coef * iqr)
  return(ret)
}

na_ratio <- function(X, by = by){
  by <- tolower(by)
  opts <- c('row', 'col', 'all')
  if (!(by %in% opts))
    stop('`by` should be one of: ', str_c(opts, collapse = ', '))
  if (by == 'all')
    ret <- sum(is.na(X))/nrow(X)/ncol(X)
  else if (by == 'row')
    ret <- apply(X, 1, function(x) sum(is.na(x)) / length(x))
  else if (by == 'col')
    ret <- apply(X, 2, function(x) sum(is.na(x)) / length(x))
  return(ret)
}
na_ratio_cutoff <- function(X, cutoff, MARGIN = 1){
  # By default, filtering rows
  na_ratio <- apply(X, MARGIN = MARGIN, function(x) sum(is.na(x)) / length(x))
  X[na_ratio < cutoff, ]
}

impute_min <- function(object, minvalue=NULL, SD=NULL, downfactor=0.8, random=T, width=0.3, seed=100) {
  #each row represents a protein;
  #each col represents a sample.
  
  if (!is.matrix(object)) object <- as.matrix(object)
  if (is.null(minvalue)) minvalue <- min(object, na.rm = T)
  nafill <- minvalue * downfactor
  if (random){
    set.seed(seed)
    if (is.null(SD)){
      x <- apply(object, 1, function(x){
        n_missing <- sum(is.na(x))
        shrinked_sd <- sd(x, na.rm = T) * width
        rnorm(n_missing, mean = nafill, sd = shrinked_sd)
      })
    } else{
      X <- data.frame(
        n_missing = apply(object, 1, function(x) sum(is.na(x)))
      )
      X$shrinked_sd <- SD * width
      x <- sapply(1:nrow(X), function(i){
        rnorm(X$n_missing[i], mean = nafill[i], sd = X$shrinked_sd[i])
      })
      # x <- apply(X, 1, function(x){
      #   rnorm(X$n_missing, mean = nafill, sd = X$shrinked_sd)
      # })
    }
    nafill <- unlist(as.vector(x))
  }
  
  object[is.na(object)] <- nafill
  return(object)
}

create_dataframe_from_list <- function(named_list, method = 'extend') {
  method_list <- c('extend', 'trim')
  if (!(method %in% method_list)){
    stop('method should be extend or trim')
  }
  if (method == 'extend'){
    max_length <- max(sapply(named_list, length))
    new_list <- lapply(named_list, function(v){
      c(v, rep(NA, max_length - length(v)))
    })
  } else if (method == 'trim'){
    min_length <- min(sapply(named_list, length))
    new_list <- lapply(named_list, function(v) head(v, min_length))
  }
  
  result_df <- as.data.frame(new_list)
  return(result_df)
}

my_color_palette <- function(n, ..., visible = T){
  require(RColorBrewer)
  col <- colorRampPalette(unlist(list(...)))(n)
  if(visible) {image(x = 1:n, y = 1, z = as.matrix(1:n), col = col)}
  return(col)
}

my_corrplot <- function(cors, color_begin = 0.2, color_end = 1, col.lim = c(0, 1.0), order = 'original'){
  require(corrplot)
  cor_colors <- viridis::viridis(101, alpha = 1, begin = color_begin, end = color_end, direction = -1, option = 'A') %>% setNames(seq(col.lim[1], col.lim[2], length.out = 101))
  corrplot(cors, order = order, type = "upper", diag = T, tl.pos = 'tl',
           cl.cex = 2, col.lim = col.lim, method = 'square', 
           tl.col = "black", addCoef.col = "#FFFFFF", number.cex = 1,
           col = cor_colors, is.corr = F) %>% print()
  corrplot(cors, order = order, type = "lower", diag = F,
           method = 'ellipse', add = T, tl.pos = 'n',
           col.lim = col.lim, col = cor_colors, is.corr = F) %>% print()
  return(NULL)
}

enrichment_analysis_env_load <- function(){
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(KEGG.db)
  library(ReactomePA)
}

# project-specific function -----------------------------------------------
pairwise_comparison <- function(df, nameCol, valueCol, test.fun, ...){
  # usage like:
  # pairwise_comparison(df, 'Type', '# protein groups', wilcox.test)
  comparisons <- combn(sort(unique(df[, nameCol])), 2, simplify = F)
  ret_ls <- lapply(comparisons, function(comp){
    x1 <- df[df[, nameCol] == comp[1], valueCol]
    x2 <- df[df[, nameCol] == comp[2], valueCol]
    tryCatch(test.fun(x1, x2, ...), error = function(e) NA)
  })
  tibble(Comparison = comparisons, test.ret = ret_ls)
}
pairwise_comparison_pvalue <- function(df, nameCol, valueCol, padj.method = 'BH', nSignif = 2, test.fun, ...){
  comparisons <- combn(sort(unique(df[, nameCol])), 2, simplify = F)
  ret_ls <- lapply(comparisons, function(comp){
    x1 <- df[df[, nameCol] == comp[1], valueCol]
    x2 <- df[df[, nameCol] == comp[2], valueCol]
    tryCatch(test.fun(x1, x2, ...), error = function(e) NA)
  })
  stat.ret <- data.frame(
    group1 = sapply(comparisons, function(comp) comp[1]),
    group2 = sapply(comparisons, function(comp) comp[2]),
    statistic = sapply(ret_ls, function(ret) ret$statistic),
    p = sapply(ret_ls, function(ret) ret$p.value)
  )
  if(!is.null(padj.method)){
    stat.ret$p.adj <- p.adjust(stat.ret$p, method = padj.method)
  }
  stat.ret[, -(1:3)] <- apply(stat.ret[, -(1:3)], 2, signif, nSignif)
  return(stat.ret)
}
create_stat_value <- function(df, nameCol, valueCol, padj.method = 'BH', nSignif = 2, test.fun, ...){
  df.stat <- df %>% 
    pairwise_comparison_pvalue(nameCol, valueCol, padj.method = padj.method, nSignif = nSignif, test.fun = test.fun, ...)
  df.stat$y.position <- max(df[, valueCol], na.rm = T) * seq(1.05, by = 0.08, length.out = nrow(df.stat))
  return(df.stat)
}

calc_pep1 <- function(pepseq){
  data.frame(Peptide = pepseq,
             pI_EMBOSS = Peptides::pI(pepseq, pKscale = 'EMBOSS'),
             pI_Lehninger = Peptides::pI(pepseq, pKscale = 'Lehninger'),
             NetCharge_EMBOSS7 = Peptides::charge(pepseq, pH = 7, pKscale = 'EMBOSS'),
             NetCharge_Lehninger7 = Peptides::charge(pepseq, pH = 7, pKscale = 'Lehninger'),
             instaIndex = Peptides::instaIndex(pepseq),
             aIndex = Peptides::aIndex(pepseq)
  )
}
calc_pep2 <- function(pepseq){
  X1 <- Peptides::crucianiProperties(pepseq) %>% setNames(pepseq) %>% as.data.frame()
  X2 <- Peptides::kideraFactors(pepseq) %>% setNames(pepseq) %>% as.data.frame()
  X3 <- Peptides::zScales(pepseq) %>% setNames(pepseq) %>% as.data.frame()
  rbind(X1, X2, X3) %>% t() %>% as.data.frame() %>% rownames_to_column('Peptide')
}
calc_missCleav <- function(pepseq){
  df_missCleav <- data.frame(pepseq = pepseq)
  df_missCleav$MissedCleavage <- 0
  df_missCleav$MissedCleavage <- apply(df_missCleav, 1, function(row){
    Seq <- row['pepseq']
    # exclude the last cleavage site
    if (stringr::str_sub(Seq, -1, -1) == 'P' & stringr::str_sub(Seq, -2, -2) %in% c('K', 'R')){
      Seq <- stringr::str_sub(Seq, 1, -3)
    }else{
      Seq <- stringr::str_sub(Seq, 1, -2)
    }
    #
    Seqs <- stringr::str_sub(Seq, 1:nchar(Seq), 1:nchar(Seq))
    is_KR <- Seqs %in% c('K', 'R')
    isnt_P <- Seqs != 'P'
    flag <- c(isnt_P[2:length(isnt_P)], T)
    mis_cleav <- is_KR & flag
    return(sum(mis_cleav))
  })
  return(df_missCleav)
}
my_pvcaBatchAssess <- function (theDataMatrix, expInfo, threshold) {
  # Cite: Bushel P (2024). pvca: Principal Variance Component Analysis (PVCA). R package version 1.44.0.
  # https://doi.org/10.1002/9780470685983.ch12
  # theDataMatrix, row as probs, col as samples;
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, 
                                           scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues/eigenValuesSum
  
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  pc_n <- ifelse(my_counter_2 < 3, 3, my_counter_2)
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  on.exit(options(op))
  effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                  ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lme4::lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                         expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                        na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(lme4::VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(lme4::getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames, matrixWtProp = set_colnames(randomEffectsMatrixWtProp, effectsNames), eigenData = eigenData))
}
