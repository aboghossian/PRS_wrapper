# Script for selecting the best PRS from a set of PRS-CSx models (varied phis)
# *** Must not be a multi-ancestry PRS ***
# Author: Andrew Boghossian

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(data.table))
suppressMessages(library(pscl))
suppressMessages(library(argparse))

# define command line args
parser <- ArgumentParser()
parser$add_argument("-p", "--pheno",
                    help = "Path to file with relevant phenotype.")
parser$add_argument("-c", "--pheno_col",
                    help = "Name of phenotype column.")
parser$add_argument("-v", "--val_prop", default = 0.5,
                    help = "Fraction of sample to use for parameter tuning.")
parser$add_argument("-d", "--in_dir", default = ".",
                    help = "Path to directory with .sscore files.")
parser$add_argument("-o", "--out_dir", default = ".",
                    help = "Path to directory to output results into.")
parser$add_argument("-n", "--name",
                    help = "Output file name.")
args <- parser$parse_args()

# phenotypes (must have FID, IID)
phenotype <- fread(args$pheno) %>%
  dplyr::select(all_of(c("FID", "IID", args$pheno_col))) %>%
  dplyr::filter(!is.na(.data[[args$pheno_col]]))

# covariates (requires Ma_Educ, MA_AGE, SEX, BOY, and PC1-10)
covars <- fread("/wynton/group/weiss/AndrewB/IMPaCT_final/IMPaCT_pheno.txt")

# model formulas from vars
model_formula <- as.formula(paste(args$pheno_col,
                                  "~ SCORE1_SUM + Ma_Educ + MA_AGE + SEX + as.factor(BOY) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
null_formula <- as.formula(paste(args$pheno_col,
                                  "~ Ma_Educ + MA_AGE + SEX + as.factor(BOY) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))

# scores
score_files <- list.files(args$in_dir, ".sscore", full.names = T)
print(score_files)

# FIDs and IIDs present in scores
individuals <- fread(score_files[[1]]) %>%
  dplyr::rename(FID = "#FID") %>%
  dplyr::select(FID, IID)

# make validation subset
val_subset <- phenotype %>%
  dplyr::sample_frac(size = as.numeric(args$val_prop)) %>%
  dplyr::inner_join(individuals)

# loop through and calculate in sample R2
scores <- list()  # stores results
print("Finding best phi value...")
for (i in 1:length(score_files)) {
  
  # join with phenotype and covariates
  score_df <- fread(score_files[[i]]) %>%
    dplyr::rename(FID = "#FID") %>%
    dplyr::inner_join(phenotype) %>%
    dplyr::inner_join(covars)
  
  # get validation set
  score_val <- score_df %>%
    dplyr::inner_join(val_subset)
  
  # within val fit
  score_fit <- glm(data = score_val,
                   formula = model_formula)
  
  # add R2 and score df to results list
  r2 <- pscl::pR2(score_fit)[["McFadden"]]  # pseudo R2 for logistic regression
  scores[[i]] <- list(score = score_df, r2 = r2)
}

# find the maximum R2 and associated score
best <- which.max(purrr::map(scores, 2) %>% unlist())
best_score <- scores[[best]]$score

# write score file to output directory with given name
write_delim(best_score %>% dplyr::select(FID, IID, SCORE1_SUM),
            file.path(args$out_dir, args$name),
            delim = "\t")
# output best R2 to console
print(paste("The best in sample R2 was",
            scores[[best]]$r2,
            "from",
            score_files[[best]]))

# fit full sample model with best PRS (phi optimized)
print("Fitting full model...")
best_full_fit <- glm(data = best_score,
                     formula = model_formula)
summary(best_full_fit)

# null fit
null_full_fit <- glm(data = best_score,
                     formula = null_formula)

# gained R2 = fill - null (what does PRS add)
print("Gained R2 in full data set was:")
pscl::pR2(best_full_fit) - pscl::pR2(null_full_fit)
