# Clear environment by unloading all loaded libraries
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Load necessary libraries
library(magrittr)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(corrplot)
library(pheatmap)
library(ggsci)
library(viridis)
# library(Peptides)
# library(lme4)

# Universal Values --------------------------------------------------------
nucleus_type_order <- c('1nucleus', 'cell_1nucleus', 'cell_2nucleus', 'blank_control')
region_order <- c('N', 'L', 'H', 'C')

# Defining color palettes for specific uses
project_colors <- c(ProteomEx = '#00008B', FAXP = '#C32022')
gel_tip_colors <- c(InGel = '#00008B', InTip = '#C32022')

cor_colors <- viridis(101, alpha = 1, begin = 0.2, end = 1, direction = -1, option = 'A') %>%
  setNames(seq(0, 1, length.out = 101))

detect_range_colors <- c('In-tip DDA-MS' = '#DF8F44FF', 
                         'In-gel DDA-MS' =  '#008280FF', 
                         'In-tip DIA-MS' = '#B24745FF', 
                         'In-gel DIA-MS' = '#00A1D5FF')

# Separate color palettes for different plots
type_colors <- setNames(c("#358DB9FF", "#CF4E9CFF", "#8C57A2FF", "#2E2A2BFF"), nucleus_type_order)
type_colors <- c(type_colors, 
                 'Intra.Group.Nucleus' = type_colors['1nucleus'],
                 'Intra.Group.Mono' = type_colors['cell_1nucleus'],
                 'Intra.Group.Bi' = type_colors['cell_2nucleus'],
                 'Inter.Groups' = type_colors['blank_control'],
                 'Nucleus' = type_colors['1nucleus'],
                 'Cell (Mononuclear)' = type_colors['cell_1nucleus'],
                 'Cell (Binuclear)' = type_colors['cell_2nucleus'])

heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(50)

region_colors <- c(N = '#000080', L = '#FFD700', H = '#CC5500', C = '#C71585')

patient_colors <- setNames(c('#77B5FE', '#4682B4', '#004488'), paste0('P', 1:3))

pca_colors <- setNames(c('#D7DA86', '#330A5FFF', '#FCB519FF', '#ED6925FF', 
                         '#781C6DFF', '#BB3754FF', '#000004FF'),
                       c('resid', 'Region', 'Patient:Region', 'Patient', 'Slide:Region', 
                         'Patient:Slide', 'Slide'))

# Universal Functions -----------------------------------------------------
cv <- function(x, na.rm = TRUE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}

get_outliers <- function(vec, coef = 1.5) {
  # outliers based on Q1 and Q3
  stats <- quantile(vec, na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  c(stats[2] - coef * iqr, stats[4] + coef * iqr)
}

na_ratio <- function(X, by = c('row', 'col', 'all')) {
  by <- match.arg(by)
  if (by == 'all') {
    sum(is.na(X)) / (nrow(X) * ncol(X))
  } else if (by == 'row') {
    rowSums(is.na(X)) / ncol(X)
  } else {
    colSums(is.na(X)) / nrow(X)
  }
}

na_ratio_cutoff <- function(X, cutoff, MARGIN = 1) {
  na_ratio <- apply(X, MARGIN, function(x) mean(is.na(x)))
  X[na_ratio < cutoff, ]
}

impute_min <- function(object, minvalue = NULL, SD = NULL, downfactor = 0.8, random = TRUE, width = 0.3, seed = 100) {
  # each row represents a protein;
  # each col represents a sample.
  
  if (!is.matrix(object)) object <- as.matrix(object)
  if (is.null(minvalue)) minvalue <- min(object, na.rm = TRUE)
  
  nafill <- minvalue * downfactor
  if (random) {
    set.seed(seed)  # Ensure reproducibility
    
    if (is.null(SD)) {
      # If SD is not provided, calculate it from the data
      x <- apply(object, 1, function(row) {
        n_missing <- sum(is.na(row))
        shrinked_sd <- sd(row, na.rm = TRUE) * width
        if (n_missing > 0) {
          rnorm(n_missing, mean = nafill, sd = shrinked_sd)
        } else {
          NULL
        }
      })
    } else {
      # If SD is provided, use it for imputation
      X <- data.frame(
        n_missing = apply(object, 1, function(row) sum(is.na(row)))
      )
      X$shrinked_sd <- SD * width
      x <- sapply(1:nrow(X), function(i) {
        if (X$n_missing[i] > 0) {
          rnorm(X$n_missing[i], mean = nafill[i], sd = X$shrinked_sd[i])
        } else {
          NULL
        }
      })
    }
    nafill_values <- unlist(as.vector(x))
  }
  
  # Impute missing values
  object[is.na(object)] <- nafill_values
  return(object)
}

create_dataframe_from_list <- function(named_list, method = c('extend', 'trim')) {
  # Validate the method parameter
  method <- match.arg(method)
  
  # Check that input is a list
  if (!is.list(named_list)) {
    stop("Input must be a list.")
  }
  
  # Calculate the lengths of each element in the list once
  lengths <- sapply(named_list, length)
  
  # Adjust the lists based on the method chosen
  if (method == 'extend') {
    max_length <- max(lengths)
    new_list <- lapply(named_list, function(v) c(v, rep(NA, max_length - length(v))))
  } else {
    min_length <- min(lengths)
    new_list <- lapply(named_list, function(v) head(v, min_length))
  }
  
  # Convert the adjusted list to a data frame
  return(as.data.frame(new_list))
}

my_color_palette <- function(n, ..., visible = TRUE) {
  require(RColorBrewer)
  col <- colorRampPalette(unlist(list(...)))(n)
  if (visible) { image(x = 1:n, y = 1, z = as.matrix(1:n), col = col) }
  return(col)
}

my_corrplot <- function(cors, color_begin = 0.2, color_end = 1, col.lim = c(0, 1.0), order = 'original') {
  require(corrplot)
  cor_colors <- viridis(101, alpha = 1, begin = color_begin, end = color_end, direction = -1, option = 'A') %>%
    setNames(seq(col.lim[1], col.lim[2], length.out = 101))
  corrplot(cors, order = order, type = "upper", diag = TRUE, tl.pos = 'tl',
           cl.cex = 2, col.lim = col.lim, method = 'square', 
           tl.col = "black", addCoef.col = "#FFFFFF", number.cex = 1,
           col = cor_colors, is.corr = FALSE) %>% print()
  corrplot(cors, order = order, type = "lower", diag = FALSE,
           method = 'ellipse', add = TRUE, tl.pos = 'n',
           col.lim = col.lim, col = cor_colors, is.corr = FALSE) %>% print()
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

# Project-Specific Functions --------------------------------------
# Helper function for pairwise comparisons
perform_pairwise_comparison <- function(df, nameCol, valueCol, test.fun, ...) {
  comparisons <- combn(sort(unique(df[[nameCol]])), 2, simplify = FALSE)
  results <- lapply(comparisons, function(comp) {
    x1 <- df[[valueCol]][df[[nameCol]] == comp[1]]
    x2 <- df[[valueCol]][df[[nameCol]] == comp[2]]
    tryCatch({
      test.fun(x1, x2, ...)
    }, error = function(e) {
      warning(sprintf("Error in test between %s and %s: %s", comp[1], comp[2], e$message))
      return(NA)
    })
  })
  list(comparisons = comparisons, results = results)
}

# Simplified pairwise comparison
pairwise_comparison <- function(df, nameCol, valueCol, test.fun, ...) {
  results <- perform_pairwise_comparison(df, nameCol, valueCol, test.fun, ...)
  tibble(Comparison = results$comparisons, test.ret = results$results)
}

# Pairwise comparison with p-value adjustment
pairwise_comparison_pvalue <- function(df, nameCol, valueCol, padj.method = 'BH', nSignif = 2, test.fun, ...) {
  results <- perform_pairwise_comparison(df, nameCol, valueCol, test.fun, ...)
  
  stat.ret <- tibble(
    group1 = sapply(results$comparisons, function(comp) comp[1]),
    group2 = sapply(results$comparisons, function(comp) comp[2]),
    statistic = sapply(results$results, function(ret) ret$statistic),
    p = sapply(results$results, function(ret) ret$p.value)
  )
  
  # Adjust p-values if a method is provided
  if (!is.null(padj.method)) {
    stat.ret <- stat.ret %>%
      mutate(p.adj = p.adjust(p, method = padj.method))
  }
  
  # Round the values to the specified number of significant digits
  stat.ret <- stat.ret %>%
    mutate(across(statistic:p.adj, ~ signif(.x, nSignif)))
  
  return(stat.ret)
}

# Create statistical values and add y-position for plotting
create_stat_value <- function(df, nameCol, valueCol, padj.method = 'BH', nSignif = 2, test.fun, ..., y_increment = 0.08) {
  df.stat <- pairwise_comparison_pvalue(df, nameCol, valueCol, padj.method = padj.method, nSignif = nSignif, test.fun = test.fun, ...)
  
  max_value <- max(df[[valueCol]], na.rm = TRUE)
  df.stat$y.position <- max_value * seq(1.05, by = y_increment, length.out = nrow(df.stat))
  
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

calc_missCleav <- function(pepseq, cleavage_rule = c('K', 'R'), exclude_end_with = 'P') {
  require(stringr)
  df_missCleav <- data.frame(pepseq = pepseq)
  
  df_missCleav$MissedCleavage <- apply(df_missCleav, 1, function(row) {
    Seq <- row['pepseq']
    
    # Exclude the last cleavage site based on the end rule
    if (str_sub(Seq, -1, -1) == exclude_end_with && 
        str_sub(Seq, -2, -2) %in% cleavage_rule) {
      Seq <- str_sub(Seq, 1, -3)
    } else {
      Seq <- str_sub(Seq, 1, -2)
    }
    
    # Determine cleavage points
    Seqs <- str_sub(Seq, 1:nchar(Seq), 1:nchar(Seq))
    is_KR <- Seqs %in% cleavage_rule
    isnt_P <- Seqs != exclude_end_with
    flag <- c(isnt_P[-1], TRUE)
    
    # Count missed cleavages
    mis_cleav <- is_KR & flag
    return(sum(mis_cleav))
  })
  
  return(df_missCleav)
}

my_pvcaBatchAssess <- function(theDataMatrix, expInfo, threshold) {
  # Cite: Bushel P (2024). pvca: Principal Variance Component Analysis (PVCA). R package version 1.44.0.
  # https://doi.org/10.1002/9780470685983.ch12
  # theDataMatrix, row as probs, col as samples;
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  theDataMatrixCentered_transposed <- apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered <- t(theDataMatrixCentered_transposed)
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues <- eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix <- eigenData$vectors
  eigenValuesSum <- sum(eigenValues)
  percents_PCs <- eigenValues / eigenValuesSum
  
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 <- 0
  my_sum_2 <- 1
  for (i in ev_n:1) {
    my_sum_2 <- my_sum_2 - percents_PCs[i]
    if (my_sum_2 <= threshold) {
      my_counter_2 <- my_counter_2 + 1
    }
  }
  pc_n <- ifelse(my_counter_2 < 3, 3, my_counter_2)
  pc_data_matrix <- matrix(0, nrow = (expDesignRowN * pc_n), ncol = 1)
  mycounter <- 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] <- eigenVectorsMatrix[j, i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- colnames(exp_design)
  for (i in seq_along(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  on.exit(options(op))
  effects_n <- expDesignColN + choose(expDesignColN, 2) + 1
  randomEffectsMatrix <- matrix(0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in seq_along(variables)) {
    mod <- paste("(1|", variables[i], ")", sep = "")
    model.func[index] <- mod
    index <- index + 1
  }
  for (i in seq_along(variables) - 1) {
    for (j in (i + 1):length(variables)) {
      mod <- paste("(1|", variables[i], ":", variables[j], ")", sep = "")
      model.func[index] <- mod
      index <- index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y <- (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lme4::lmer(funct, Data[y:(((i - 1) * expDesignRowN) + expDesignRowN), ], REML = TRUE, verbose = FALSE, na.action = na.omit)
    randomEffectsMatrix[i, ] <- c(unlist(lme4::VarCorr(Rm1ML)), resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(lme4::getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n) {
    mySum <- sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] <- randomEffectsMatrix[i, j] / mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n) {
    weight <- eigenValues[i] / eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] <- randomEffectsMatrixStdze[i, j] * weight
    }
  }
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum <- sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(0, nrow = 1, ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] <- randomEffectsSums[j] / totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames, matrixWtProp = set_colnames(randomEffectsMatrixWtProp, effectsNames), eigenData = eigenData))
}
