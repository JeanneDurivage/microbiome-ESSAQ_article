Differential abundance analyses of 16S and ITS ASVs response to tillage
and organic fertilisation
================

- [1 Setup](#1-setup)
- [2 Function creation](#2-function-creation)
- [3 Data preparation](#3-data-preparation)
- [4 16S analysis](#4-16s-analysis)
  - [4.1 Humic orthic gleysols](#41-humic-orthic-gleysols)
  - [4.2 Orthic gleysols](#42-orthic-gleysols)
  - [4.3 Humo-ferric podzol](#43-humo-ferric-podzol)
- [5 ITS analysis](#5-its-analysis)
  - [5.1 Humic orthic gleysols](#51-humic-orthic-gleysols)
  - [5.2 Orthic gleysols](#52-orthic-gleysols)
  - [5.3 Humo-ferric podzol](#53-humo-ferric-podzol)
- [6 Results visualisation](#6-results-visualisation)
- [7 Session info](#7-session-info)

**This example codes shows the ANCOM-BC2 analyses for the 16S and the
ITS data.** Differential abundance analysis of prokaryotic C cycling
potential is done in another script.

# 1 Setup

``` r
set.seed(241) # random.org
```

This R projects uses `renv`. A `.lock` file to restore adequate package
versions is available as part of the GitHub. See the [package
vignette](https://rstudio.github.io/renv/articles/renv.html#collaboration)
to know more.

``` r
library(phyloseq)    # accepted data format for ANCOM-BC2()
library(wesanderson) # plot colors
```

    ## Warning: le package 'wesanderson' a été compilé avec la version R 4.4.3

``` r
library(ggpubr)      # assembling plots
```

    ## Le chargement a nécessité le package : ggplot2

    ## Warning: le package 'ggplot2' a été compilé avec la version R 4.4.3

``` r
library(ANCOMBC)
```

    ## Warning: le package 'ANCOMBC' a été compilé avec la version R 4.4.2

    ## Registered S3 method overwritten by 'DescTools':
    ##   method        from       
    ##   print.palette wesanderson

    ## Registered S3 method overwritten by 'lme4':
    ##   method           from
    ##   na.action.merMod car

``` r
library(doRNG)       # used by ANCOMBC2()
```

    ## Warning: le package 'doRNG' a été compilé avec la version R 4.4.3

    ## Le chargement a nécessité le package : foreach

    ## Warning: le package 'foreach' a été compilé avec la version R 4.4.3

    ## Le chargement a nécessité le package : rngtools

    ## Warning: le package 'rngtools' a été compilé avec la version R 4.4.3

``` r
library(foreach)     # used by ANCOMBC2()
library(rngtools)    # used by ANCOMBC2()
library(tidyverse) 
```

    ## Warning: le package 'tibble' a été compilé avec la version R 4.4.3

    ## Warning: le package 'tidyr' a été compilé avec la version R 4.4.3

    ## Warning: le package 'readr' a été compilé avec la version R 4.4.3

    ## Warning: le package 'purrr' a été compilé avec la version R 4.4.3

    ## Warning: le package 'dplyr' a été compilé avec la version R 4.4.3

    ## Warning: le package 'stringr' a été compilé avec la version R 4.4.3

    ## Warning: le package 'forcats' a été compilé avec la version R 4.4.3

    ## Warning: le package 'lubridate' a été compilé avec la version R 4.4.3

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ lubridate 1.9.5     ✔ tibble    3.3.1
    ## ✔ purrr     1.2.2     ✔ tidyr     1.3.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ purrr::accumulate() masks foreach::accumulate()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ purrr::when()       masks foreach::when()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
theme_set(theme_bw()) # Default ggplot theme

library(conflicted)
conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename)
```

    ## [conflicted] Will prefer dplyr::filter over any other package.
    ## [conflicted] Will prefer dplyr::select over any other package.
    ## [conflicted] Will prefer dplyr::rename over any other package.

``` r
github_path       <- "C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article" # Change me!
fig_path          <- file.path(github_path, "figures")
article_fig_path  <- file.path(github_path, "figures_article")
data_path         <- file.path(github_path, "data")
private_data_path <- file.path(github_path, "data_private")

knitr::opts_chunk$set(fig.path = paste0(fig_path, "/"))
```

# 2 Function creation

Function to extract significant taxa from the ANCOM-BC2.

``` r
extract_signif_taxa <- function(ancombc_output, community, ped_gr, tax_table) {
  # ancombc_output : ANCOM-BC2 analysis result
  # community : Either "16S" of "ITS"
  # ped_gr : Either "HOG", "OG" or "HFP"
  # tax_table : taxonomic table

  res_OF <-
    ancombc_output$res_pair |>
    select(taxon, starts_with("lfc_"), 
           starts_with("q_"), starts_with("passed_ss_")) |>
    pivot_longer(
      cols = -taxon, names_to = c(".value", "variable"),
      names_pattern = "(lfc|q|passed_ss)_(.*)"
    ) |>
    filter(q <= pval & passed_ss == TRUE) |>
    separate(variable, into = c("Comp1", "Comp2"), 
             sep = "_OF_type", fill = "right") |>
    mutate(
      signf = case_when(q <= 0.001 / 3 ~ "***", 
                        q <= 0.01 / 3 ~ "**", 
                        q <= 0.05 / 3 ~ "*"),
      Comp1 = gsub("OF_type", "", Comp1), 
      Comp2 = gsub("OF_type", "", Comp2),
      Comp2 = if_else(is.na(Comp2), "No-OF", Comp2)
    ) |>
    filter(Comp1 != "Other-OF" & Comp2 != "Other-OF") |> # No interested in these comparisons
    # adding metadata
    mutate(
      Community = community, 
      Ped._gr. = ped_gr,
      Treatment = paste0(Comp1, " vs ", Comp2)
    ) |> 
    select(-passed_ss, -Comp1, -Comp2) 

  cat("There are", length(unique(res_OF$taxon)), 
      "unique significant taxa for the paiwise OF comparisons.\n")

  res_STIR <-
    ancombc_output$res |>
    select(taxon, starts_with("lfc_"), 
           starts_with("q_"), starts_with("passed_ss_")) |>
    pivot_longer(
      cols = -taxon, names_to = c(".value", "variable"),
      names_pattern = "(lfc|q|passed_ss)_(.*)"
    ) |>
    filter(variable == "STIR") |>
    filter(q <= pval & passed_ss == TRUE) |>
    mutate(signf = case_when(q <= 0.001 / 3 ~ "***", 
                             q <= 0.01 / 3 ~ "**", 
                             q <= 0.05 / 3 ~ "*")) |>
    select(-passed_ss) |>
    rename(Treatment = variable) |> 
    # adding metadata
    mutate(Community = community, 
           Ped._gr. = ped_gr) 


  cat("There are", length(unique(res_STIR$taxon)), 
      "unique significant taxa for STIR.")

  # Putting everything together
  res_all <- rbind(res_OF, res_STIR) |>
    # Adding the phylum
    left_join(tax_table |> select(Phylum, Family) |> unique(), 
              by = join_by(taxon == Family)) |>
    mutate(
      taxon = str_remove(taxon, "f__"),
      Phylum = str_remove(Phylum, "p__"),
      taxon_phylum = paste0(taxon, " (", Phylum, ")"),
      fold_change = if_else(lfc < 0, -exp(abs(lfc)), exp(lfc))
    ) |>
    select(-Phylum) 

  return(res_all)
}
```

# 3 Data preparation

For confidientality reasons, the raw data is not publically available.

``` r
# Importation
otu_16S <- read.csv(file.path(private_data_path, "otu_16S.csv"), row.names = 1)
otu_ITS <- read.csv(file.path(private_data_path, "otu_ITS.csv"), row.names = 1)
tax_16S <- read.csv(file.path(private_data_path, "tax_16S.csv"), row.names = 1)
tax_ITS <- read.csv(file.path(private_data_path, "tax_ITS.csv"), row.names = 1)
env_raw <- read.csv(file.path(private_data_path, "environment_data.csv"), row.names = 1)
```

This scripts uses 5 datasets :

- `env_raw`, which contains all necessary soil-, climate- and
  farming-related variables.

- `otu_16S` and `otu_ITS`, which contains the ASV abundance table from
  the 16S and ITS rRNA amplicon sequencing

- `tax_16S` and `tax_ITS`, which contains the taxonomic attribution to
  each ASV.

We will prepare the `env` table by selecting necessary variables and
converting character vectors into ordered factors.

``` r
env_ASV <- env_raw |> select(OF_type, STIR, Ped._gr.) |> 
                      mutate(
                             # Scaling of numeric response variables
                             STIR = as.numeric(scale(STIR)),
                             # Logigal ordering of factor levels
                             Ped._gr. = factor(Ped._gr., levels = c("Humo-ferric podzol",
                                                                    "Orthic gleysol",
                                                                    "Humic orthic gleysol")),
                             OF_type  = factor(OF_type, levels = c("No-OF",
                                                                   "L-Hog",
                                                                   "L-Dairy",
                                                                   "S-Poultry",
                                                                   "Other-OF"))
                             ) 
```

Overview of the environmental data :

``` r
skimr::skim(env_ASV)
```

|                                                  |         |
|:-------------------------------------------------|:--------|
| Name                                             | env_ASV |
| Number of rows                                   | 80      |
| Number of columns                                | 3       |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |         |
| Column type frequency:                           |         |
| factor                                           | 2       |
| numeric                                          | 1       |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |         |
| Group variables                                  | None    |

Data summary

**Variable type: factor**

| skim_variable | n_missing | complete_rate | ordered | n_unique | top_counts |
|:---|---:|---:|:---|---:|:---|
| OF_type | 0 | 1 | FALSE | 5 | No-: 22, Oth: 22, L-D: 17, L-H: 12 |
| Ped.\_gr. | 0 | 1 | FALSE | 3 | Hum: 46, Hum: 21, Ort: 13 |

**Variable type: numeric**

| skim_variable | n_missing | complete_rate | mean |  sd |    p0 |   p25 |  p50 |  p75 | p100 | hist  |
|:--------------|----------:|--------------:|-----:|----:|------:|------:|-----:|-----:|-----:|:------|
| STIR          |         0 |             1 |    0 |   1 | -1.51 | -0.82 | 0.04 | 0.71 | 4.17 | ▇▇▅▁▁ |

- `OF_type` corresponds to the type of organic fertilizer.

- `STIR` is the tillage intensity index.

- `Ped._gr.` is the pedologic group.

``` r
table(env_ASV$OF_type, env_ASV$Ped._gr.)
```

    ##            
    ##             Humo-ferric podzol Orthic gleysol Humic orthic gleysol
    ##   No-OF                      3              3                   16
    ##   L-Hog                      4              2                    6
    ##   L-Dairy                    7              4                    6
    ##   S-Poultry                  0              3                    4
    ##   Other-OF                   7              1                   14

The `Other-OF` category will be a problem in `Orthic gleysol` where its
sample size `= 1`.

This sample will be removed.

``` r
env_ASV <- env_ASV |> filter(!OF_type == "Other-OF" | !Ped._gr. == "Orthic gleysol")
```

Creation of `phyloseq` objects :

``` r
phy_16S <- phyloseq(otu_table(otu_16S, taxa_are_rows = TRUE),
                    tax_table(as.matrix(tax_16S)),
                    sample_data(env_ASV))
phy_ITS <- phyloseq(otu_table(otu_ITS, taxa_are_rows = TRUE),
                    tax_table(as.matrix(tax_ITS)),
                    sample_data(env_ASV))
```

# 4 16S analysis

Determining the proportion of ASVs with a taxonomic attribution.

``` r
n_OTU_att_16S <- phy_16S@tax_table |> 
  as.data.frame() |> 
  mutate(Family = ifelse(Family == "f__", NA, Family)) |> 
  select(Family) |> 
  drop_na() |> 
  nrow() 
  
prc_OTU_att_16S <- round((n_OTU_att_16S/nrow(phy_16S@tax_table)*100), 1)

cat(prc_OTU_att_16S, "% of 16S ASVs have a taxonomic attribution at the family level.")
```

    ## 84.5 % of 16S ASVs have a taxonomic attribution at the family level.

Analyses are run by pedological group, individually. We will split the
datasets according to the pedologic groups.

``` r
grps_16S <- levels(sample_data(phy_16S)$Ped._gr.)

# Splitting 16S
phy_list_16S <- setNames(
  lapply(grps_16S, function(g) prune_samples(sample_data(phy_16S)$Ped._gr. == g, phy_16S)),
  grps_16S)
```

Setting a lower p-value to account for repeating the ANCOM-BC2 three
times on 16S data.

``` r
pval <- 0.05 / 3 
```

## 4.1 Humic orthic gleysols

``` r
output_bac_HOG <- ancombc2(data = phy_list_16S[["Humic orthic gleysol"]], 
                           tax_level = "Family",
                           fix_formula = "STIR + OF_type", 
                           group = "OF_type", pairwise = TRUE,
                           prv_cut = 5/ncol(phy_list_16S[["Humic orthic gleysol"]]@otu_table),
                           struc_zero = TRUE
                           )
```

    ## Checking the input data type ...

    ## The input data is of type: phyloseq

    ## PASS

    ## Checking the sample metadata ...

    ## The specified variables in the formula: STIR, OF_type

    ## The available variables in the sample metadata: OF_type, STIR, Ped._gr.

    ## PASS

    ## Checking other arguments ...

    ## The number of groups of interest is: 5

    ## The sample size per group is: No-OF = 16, L-Hog = 6, L-Dairy = 6, S-Poultry = 4, Other-OF = 14

    ## Warning: Small sample size detected for the following group(s): 
    ## S-Poultry
    ## Variance estimation would be unstable when the sample size is < 5 per group

    ## PASS

    ## Obtaining initial estimates ...

    ## Estimating sample-specific biases ...

    ## Conducting sensitivity analysis for pseudo-count addition to 0s ...
    ## For taxa that are significant but do not pass the sensitivity analysis,
    ## please flag them and proceed with caution, as they are likely false positives.
    ## For detailed instructions on performing sensitivity analysis,
    ## please refer to the package vignette.

    ## ANCOM-BC2 primary results ...

    ## ANCOM-BC2 multiple pairwise comparisons ...

Extracting the significant taxa.

``` r
bac_HOG <- extract_signif_taxa(output_bac_HOG, "16S", 
                               "Humic orthic gleysol", tax_16S)
```

    ## There are 6 unique significant taxa for the paiwise OF comparisons.
    ## There are 0 unique significant taxa for STIR.

## 4.2 Orthic gleysols

``` r
output_bac_OG <- ancombc2(data = phy_list_16S[["Orthic gleysol"]], 
                           tax_level = "Family",
                           fix_formula = "STIR + OF_type", 
                           group = "OF_type", pairwise = TRUE,
                           prv_cut = 5/ncol(phy_list_16S[["Orthic gleysol"]]@otu_table),
                           struc_zero = TRUE
                           )
```

    ## Checking the input data type ...

    ## The input data is of type: phyloseq

    ## PASS

    ## Checking the sample metadata ...

    ## The specified variables in the formula: STIR, OF_type

    ## The available variables in the sample metadata: OF_type, STIR, Ped._gr.

    ## PASS

    ## Checking other arguments ...

    ## The number of groups of interest is: 4

    ## The sample size per group is: No-OF = 3, L-Hog = 2, L-Dairy = 4, S-Poultry = 3

    ## Warning: Small sample size detected for the following group(s): 
    ## No-OF, L-Hog, L-Dairy, S-Poultry
    ## Variance estimation would be unstable when the sample size is < 5 per group

    ## PASS

    ## Obtaining initial estimates ...

    ## Estimating sample-specific biases ...

    ## Conducting sensitivity analysis for pseudo-count addition to 0s ...
    ## For taxa that are significant but do not pass the sensitivity analysis,
    ## please flag them and proceed with caution, as they are likely false positives.
    ## For detailed instructions on performing sensitivity analysis,
    ## please refer to the package vignette.

    ## ANCOM-BC2 primary results ...

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

    ## ANCOM-BC2 multiple pairwise comparisons ...

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

Extracting the significant taxa.

``` r
bac_OG <- extract_signif_taxa(output_bac_OG, "16S", 
                              "Orthic gleysol", tax_16S)
```

    ## There are 13 unique significant taxa for the paiwise OF comparisons.
    ## There are 0 unique significant taxa for STIR.

## 4.3 Humo-ferric podzol

``` r
output_bac_HFP <- ancombc2(data = phy_list_16S[["Humo-ferric podzol"]], 
                           tax_level = "Family",
                           fix_formula = "STIR + OF_type", 
                           group = "OF_type", pairwise = TRUE,
                           prv_cut = 5/ncol(phy_list_16S[["Humo-ferric podzol"]]@otu_table),
                           struc_zero = TRUE
                           )
```

    ## Checking the input data type ...

    ## The input data is of type: phyloseq

    ## PASS

    ## Checking the sample metadata ...

    ## The specified variables in the formula: STIR, OF_type

    ## The available variables in the sample metadata: OF_type, STIR, Ped._gr.

    ## PASS

    ## Checking other arguments ...

    ## The number of groups of interest is: 4

    ## The sample size per group is: No-OF = 3, L-Hog = 4, L-Dairy = 7, Other-OF = 7

    ## Warning: Small sample size detected for the following group(s): 
    ## No-OF, L-Hog
    ## Variance estimation would be unstable when the sample size is < 5 per group

    ## PASS

    ## Obtaining initial estimates ...

    ## Estimating sample-specific biases ...

    ## Conducting sensitivity analysis for pseudo-count addition to 0s ...
    ## For taxa that are significant but do not pass the sensitivity analysis,
    ## please flag them and proceed with caution, as they are likely false positives.
    ## For detailed instructions on performing sensitivity analysis,
    ## please refer to the package vignette.

    ## ANCOM-BC2 primary results ...

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

    ## ANCOM-BC2 multiple pairwise comparisons ...

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

Extracting the significant taxa.

``` r
bac_HFP <- extract_signif_taxa(output_bac_HFP, "16S", 
                               "Humo-ferric podzol", tax_16S)
```

    ## There are 0 unique significant taxa for the paiwise OF comparisons.
    ## There are 0 unique significant taxa for STIR.

# 5 ITS analysis

Determining the proportion of ASVs with a taxonomic attribution.

``` r
n_OTU_att_ITS <- phy_ITS@tax_table |> 
  as.data.frame() |> 
  mutate(Family = ifelse(Family == "f__", NA, Family)) |> 
  select(Family) |> 
  drop_na() |> 
  nrow() 
  
prc_OTU_att_ITS <- round((n_OTU_att_ITS/nrow(phy_ITS@tax_table)*100), 1)

cat(prc_OTU_att_ITS, "% of ITS ASVs have a taxonomic attribution at the family level.")
```

    ## 28.1 % of ITS ASVs have a taxonomic attribution at the family level.

Analyses are run by pedological group, individually. We will split the
datasets according to the pedologic groups.

``` r
grps_ITS <- levels(sample_data(phy_ITS)$Ped._gr.)

# Splitting ITS
phy_list_ITS <- setNames(
  lapply(grps_ITS, function(g) prune_samples(sample_data(phy_ITS)$Ped._gr. == g, phy_ITS)),
  grps_ITS)
```

Setting a lower p-value to account for repeating the ANCOM-BC2 three
times on ITS data.

``` r
pval <- 0.05 / 3 
```

## 5.1 Humic orthic gleysols

``` r
output_fung_HOG <- ancombc2(data = phy_list_ITS[["Humic orthic gleysol"]], 
                           tax_level = "Family",
                           fix_formula = "STIR + OF_type", 
                           group = "OF_type", pairwise = TRUE,
                           prv_cut = 5/ncol(phy_list_ITS[["Humic orthic gleysol"]]@otu_table),
                           struc_zero = TRUE
                           )
```

    ## Checking the input data type ...

    ## The input data is of type: phyloseq

    ## PASS

    ## Checking the sample metadata ...

    ## The specified variables in the formula: STIR, OF_type

    ## The available variables in the sample metadata: OF_type, STIR, Ped._gr.

    ## PASS

    ## Checking other arguments ...

    ## The number of groups of interest is: 5

    ## The sample size per group is: No-OF = 16, L-Hog = 6, L-Dairy = 6, S-Poultry = 4, Other-OF = 14

    ## Warning: Small sample size detected for the following group(s): 
    ## S-Poultry
    ## Variance estimation would be unstable when the sample size is < 5 per group

    ## PASS

    ## Obtaining initial estimates ...

    ## Estimating sample-specific biases ...

    ## Conducting sensitivity analysis for pseudo-count addition to 0s ...
    ## For taxa that are significant but do not pass the sensitivity analysis,
    ## please flag them and proceed with caution, as they are likely false positives.
    ## For detailed instructions on performing sensitivity analysis,
    ## please refer to the package vignette.

    ## ANCOM-BC2 primary results ...

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

    ## ANCOM-BC2 multiple pairwise comparisons ...

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

Extracting the significant taxa.

``` r
fung_HOG <- extract_signif_taxa(output_fung_HOG, "ITS", 
                               "Humic orthic gleysol", tax_ITS)
```

    ## There are 3 unique significant taxa for the paiwise OF comparisons.
    ## There are 0 unique significant taxa for STIR.

## 5.2 Orthic gleysols

``` r
output_fung_OG <- ancombc2(data = phy_list_ITS[["Orthic gleysol"]], 
                           tax_level = "Family",
                           fix_formula = "STIR + OF_type", 
                           group = "OF_type", pairwise = TRUE,
                           prv_cut = 5/ncol(phy_list_ITS[["Orthic gleysol"]]@otu_table),
                           struc_zero = TRUE
                           )
```

    ## Checking the input data type ...

    ## The input data is of type: phyloseq

    ## PASS

    ## Checking the sample metadata ...

    ## The specified variables in the formula: STIR, OF_type

    ## The available variables in the sample metadata: OF_type, STIR, Ped._gr.

    ## PASS

    ## Checking other arguments ...

    ## The number of groups of interest is: 4

    ## The sample size per group is: No-OF = 3, L-Hog = 2, L-Dairy = 4, S-Poultry = 3

    ## Warning: Small sample size detected for the following group(s): 
    ## No-OF, L-Hog, L-Dairy, S-Poultry
    ## Variance estimation would be unstable when the sample size is < 5 per group

    ## PASS

    ## Obtaining initial estimates ...

    ## Estimating sample-specific biases ...

    ## Conducting sensitivity analysis for pseudo-count addition to 0s ...
    ## For taxa that are significant but do not pass the sensitivity analysis,
    ## please flag them and proceed with caution, as they are likely false positives.
    ## For detailed instructions on performing sensitivity analysis,
    ## please refer to the package vignette.

    ## ANCOM-BC2 primary results ...

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

    ## ANCOM-BC2 multiple pairwise comparisons ...

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

Extracting the significant taxa.

``` r
fung_OG <- extract_signif_taxa(output_fung_OG, "ITS", 
                              "Orthic gleysol", tax_ITS)
```

    ## There are 3 unique significant taxa for the paiwise OF comparisons.
    ## There are 0 unique significant taxa for STIR.

## 5.3 Humo-ferric podzol

``` r
output_fung_HFP <- ancombc2(data = phy_list_ITS[["Humo-ferric podzol"]], 
                           tax_level = "Family",
                           fix_formula = "STIR + OF_type", 
                           group = "OF_type", pairwise = TRUE,
                           prv_cut = 5/ncol(phy_list_ITS[["Humo-ferric podzol"]]@otu_table),
                           struc_zero = TRUE
                           )
```

    ## Checking the input data type ...

    ## The input data is of type: phyloseq

    ## PASS

    ## Checking the sample metadata ...

    ## The specified variables in the formula: STIR, OF_type

    ## The available variables in the sample metadata: OF_type, STIR, Ped._gr.

    ## PASS

    ## Checking other arguments ...

    ## The number of groups of interest is: 4

    ## The sample size per group is: No-OF = 3, L-Hog = 4, L-Dairy = 7, Other-OF = 7

    ## Warning: Small sample size detected for the following group(s): 
    ## No-OF, L-Hog
    ## Variance estimation would be unstable when the sample size is < 5 per group

    ## PASS

    ## Obtaining initial estimates ...

    ## Estimating sample-specific biases ...

    ## Conducting sensitivity analysis for pseudo-count addition to 0s ...
    ## For taxa that are significant but do not pass the sensitivity analysis,
    ## please flag them and proceed with caution, as they are likely false positives.
    ## For detailed instructions on performing sensitivity analysis,
    ## please refer to the package vignette.

    ## ANCOM-BC2 primary results ...

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

    ## ANCOM-BC2 multiple pairwise comparisons ...

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## TRUE): Production de NaN

    ## Warning in pf(W_global, df1 = length(beta_hat_sub_i), df2 = dof_i, lower.tail =
    ## FALSE): Production de NaN

    ## Warning in pt(abs(W), df = dof, lower.tail = FALSE): Production de NaN

Extracting the significant taxa.

``` r
fung_HFP <- extract_signif_taxa(output_fung_HFP, "ITS", 
                               "Humo-ferric podzol", tax_ITS)
```

    ## There are 1 unique significant taxa for the paiwise OF comparisons.
    ## There are 2 unique significant taxa for STIR.

# 6 Results visualisation

``` r
df_list <- list(bac_HOG,  bac_OG,  bac_HFP,
                fung_HOG, fung_OG, fung_HFP)

df_list_not_empty <- df_list |> keep(~ nrow(.) > 0)

recap_df <- reduce(df_list_not_empty, full_join)
```

    ## Joining with `by = join_by(taxon, lfc, q, signf, Community, Ped._gr.,
    ## Treatment, taxon_phylum, fold_change)`
    ## Joining with `by = join_by(taxon, lfc, q, signf, Community, Ped._gr.,
    ## Treatment, taxon_phylum, fold_change)`
    ## Joining with `by = join_by(taxon, lfc, q, signf, Community, Ped._gr.,
    ## Treatment, taxon_phylum, fold_change)`
    ## Joining with `by = join_by(taxon, lfc, q, signf, Community, Ped._gr.,
    ## Treatment, taxon_phylum, fold_change)`

``` r
recap_df$signf <- factor(recap_df$signf, levels = c("*", "**", "***"))

all_signif_taxa <- 
  recap_df |> 
  mutate(
    Effect_size_log = signif(lfc, 3),
    Effect_size     = signif(fold_change, 3),
    # correcting the name of the taxa assigned at the order or class level
    taxon_phylum    = case_when(
      taxon_phylum == "k__Fungi_p__Basidiomycota_c__Agaricomycetes_o__Cantharellales_ (NA)" ~ 
        "Cantharellales order (Basidiomycota)", 
      taxon_phylum == "k__Fungi_p__Ascomycota_c__Eurotiomycetes_o__Chaetothyriales_ (NA)" ~ 
        "Chaetothyriales order (Ascomycota)", 
      taxon_phylum == "d__Bacteria_p__Actinobacteriota_c__Thermoleophilia_o___ (NA)" ~
        "Thermoleophilia class (Actinobacteriota)",
      .default = taxon_phylum
      )
  ) |> 
  select(Community, Ped._gr., taxon_phylum, Treatment, 
         Effect_size_log, Effect_size, signf) |> 
  arrange(Community, Ped._gr., Treatment)
```

``` r
knitr::kable(all_signif_taxa)
```

| Community | Ped.\_gr. | taxon_phylum | Treatment | Effect_size_log | Effect_size | signf |
|:---|:---|:---|:---|---:|---:|:---|
| 16S | Humic orthic gleysol | Thermoleophilia class (Actinobacteriota) | L-Dairy vs No-OF | -2.00 | -7.36 | \*\* |
| 16S | Humic orthic gleysol | Microscillaceae (Bacteroidota) | L-Hog vs No-OF | -1.27 | -3.57 | \* |
| 16S | Humic orthic gleysol | Bdellovibrionaceae (Bdellovibrionota) | L-Hog vs No-OF | -1.68 | -5.38 | \* |
| 16S | Humic orthic gleysol | Cellvibrionaceae (Proteobacteria) | L-Hog vs No-OF | -2.26 | -9.58 | \* |
| 16S | Humic orthic gleysol | Peptostreptococcaceae (Firmicutes) | S-Poultry vs L-Hog | -2.77 | -16.00 | \*\*\* |
| 16S | Humic orthic gleysol | Flavobacteriaceae (Bacteroidota) | S-Poultry vs No-OF | -1.18 | -3.26 | \* |
| 16S | Orthic gleysol | Micromonosporaceae (Actinobacteriota) | L-Dairy vs No-OF | -2.66 | -14.40 | \* |
| 16S | Orthic gleysol | Paenibacillaceae (Firmicutes) | L-Dairy vs No-OF | 2.78 | 16.10 | \* |
| 16S | Orthic gleysol | WD2101_soil_group (Planctomycetota) | L-Dairy vs No-OF | -2.15 | -8.55 | \* |
| 16S | Orthic gleysol | Azospirillaceae (Proteobacteria) | L-Dairy vs No-OF | -2.56 | -12.90 | \* |
| 16S | Orthic gleysol | A0839 (Proteobacteria) | L-Dairy vs No-OF | 4.41 | 82.40 | \* |
| 16S | Orthic gleysol | Micromonosporaceae (Actinobacteriota) | L-Hog vs No-OF | -2.46 | -11.70 | \* |
| 16S | Orthic gleysol | Roseiflexaceae (Chloroflexi) | L-Hog vs No-OF | -2.16 | -8.66 | \*\* |
| 16S | Orthic gleysol | Paenibacillaceae (Firmicutes) | L-Hog vs No-OF | 2.53 | 12.60 | \*\*\* |
| 16S | Orthic gleysol | WD2101_soil_group (Planctomycetota) | L-Hog vs No-OF | -2.22 | -9.24 | \* |
| 16S | Orthic gleysol | Hyphomonadaceae (Proteobacteria) | L-Hog vs No-OF | 2.17 | 8.75 | \* |
| 16S | Orthic gleysol | A0839 (Proteobacteria) | L-Hog vs No-OF | 4.40 | 81.70 | \*\* |
| 16S | Orthic gleysol | Xanthomonadaceae (Proteobacteria) | L-Hog vs No-OF | 1.07 | 2.91 | \* |
| 16S | Orthic gleysol | Pedosphaeraceae (Verrucomicrobiota) | L-Hog vs No-OF | -1.38 | -3.97 | \* |
| 16S | Orthic gleysol | Micromonosporaceae (Actinobacteriota) | S-Poultry vs L-Dairy | 2.60 | 13.40 | \* |
| 16S | Orthic gleysol | RBG-13-54-9 (Chloroflexi) | S-Poultry vs L-Dairy | -2.39 | -10.90 | \* |
| 16S | Orthic gleysol | WX65 (Methylomirabilota) | S-Poultry vs L-Dairy | -2.46 | -11.70 | \* |
| 16S | Orthic gleysol | WD2101_soil_group (Planctomycetota) | S-Poultry vs L-Dairy | 2.26 | 9.58 | \* |
| 16S | Orthic gleysol | Azospirillaceae (Proteobacteria) | S-Poultry vs L-Dairy | 2.71 | 15.00 | \* |
| 16S | Orthic gleysol | A0839 (Proteobacteria) | S-Poultry vs L-Dairy | -3.70 | -40.60 | \*\* |
| 16S | Orthic gleysol | Micromonosporaceae (Actinobacteriota) | S-Poultry vs L-Hog | 2.39 | 10.90 | \* |
| 16S | Orthic gleysol | Roseiflexaceae (Chloroflexi) | S-Poultry vs L-Hog | 2.46 | 11.80 | \*\* |
| 16S | Orthic gleysol | WD2101_soil_group (Planctomycetota) | S-Poultry vs L-Hog | 2.34 | 10.40 | \*\* |
| 16S | Orthic gleysol | Isosphaeraceae (Planctomycetota) | S-Poultry vs L-Hog | 2.95 | 19.10 | \*\* |
| 16S | Orthic gleysol | A0839 (Proteobacteria) | S-Poultry vs L-Hog | -3.70 | -40.30 | \*\* |
| 16S | Orthic gleysol | Streptosporangiaceae (Actinobacteriota) | S-Poultry vs No-OF | 1.98 | 7.28 | \* |
| ITS | Humic orthic gleysol | Mycosphaerellaceae (Ascomycota) | L-Dairy vs No-OF | -2.58 | -13.20 | \* |
| ITS | Humic orthic gleysol | Cantharellales order (Basidiomycota) | L-Dairy vs No-OF | -2.57 | -13.10 | \* |
| ITS | Humic orthic gleysol | Sporormiaceae (Ascomycota) | S-Poultry vs L-Dairy | -2.54 | -12.70 | \* |
| ITS | Humic orthic gleysol | Cantharellales order (Basidiomycota) | S-Poultry vs L-Dairy | 3.91 | 49.90 | \*\* |
| ITS | Humo-ferric podzol | Strophariaceae (Basidiomycota) | L-Dairy vs No-OF | -4.22 | -68.00 | \*\*\* |
| ITS | Humo-ferric podzol | Chaetothyriales order (Ascomycota) | STIR | -1.19 | -3.27 | \* |
| ITS | Humo-ferric podzol | Cephalothecaceae (Ascomycota) | STIR | 1.84 | 6.29 | \* |
| ITS | Orthic gleysol | Sanchytriaceae (Monoblepharomycota) | L-Dairy vs No-OF | 4.38 | 79.90 | \* |
| ITS | Orthic gleysol | Sanchytriaceae (Monoblepharomycota) | S-Poultry vs L-Dairy | -5.06 | -157.00 | \* |
| ITS | Orthic gleysol | Leotiaceae (Ascomycota) | S-Poultry vs L-Hog | 4.76 | 116.00 | \* |
| ITS | Orthic gleysol | Onygenaceae (Ascomycota) | S-Poultry vs No-OF | -4.33 | -76.00 | \*\*\* |

Saving the results (corresponds to suplementary material 8)

``` r
write.csv2(all_signif_taxa, file = file.path(data_path, "ANCOMBC.csv"), row.names = FALSE)
```

Plotting results

``` r
# Graph data, available on the GitHub
all_signif_taxa       <- read.csv2(file.path(data_path, "ANCOMBC.csv"))
```

``` r
p_16S <- 
  all_signif_taxa |>
  filter(Community == "16S") |>
  # Adding empty rows for pedologic groups and treatments without significant taxa
  add_row(
    Community       = rep("16S", 3),
    Ped._gr.        = rep("Humo-ferric podzol", 3),
    taxon_phylum    = rep(" ", 3),
    Treatment       = c("STIR", "S-Poultry vs No-OF", "L-Dairy vs L-Hog"),
    Effect_size_log = rep(0, 3),
    Effect_size     = rep(0, 3),
    signf           = rep("***", 3)
  ) |>
  ggplot(aes(x = taxon_phylum, y = Effect_size_log, fill = signf)) +
  facet_grid(Ped._gr. ~ Treatment, scales = "free_x") +
  geom_col(show.legend = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "p-value", y = "Log-fold change", x = "Taxonomic family") +
  scale_fill_manual(values = c("#95D840", "#33638D", "#481567"), 
                    drop = FALSE,
                    labels=c('<= 0.017', '<= 0.0033', '<= 0.00033')) +
  geom_hline(yintercept = 0, colour = "firebrick4") +
  theme(plot.margin = unit(c(1,1,1,2), "cm"))
```

``` r
p_ITS <- 
  all_signif_taxa |>
  filter(Community == "ITS") |>
  # Adding empty rows for pedologic groups and treatments without significant taxa
  add_row(
    Community       = rep("ITS", 3),
    Ped._gr.        = rep("Humo-ferric podzol", 3),
    taxon_phylum    = rep(" ", 3),
    Treatment       = c("L-Hog vs No-OF", "S-Poultry vs L-Hog", "L-Dairy vs L-Hog"),
    Effect_size_log = rep(0, 3),
    Effect_size     = rep(0, 3),
    signf           = rep("***", 3)
  ) |>
  ggplot(aes(x = taxon_phylum, y = Effect_size_log, fill = signf)) +
  facet_grid(Ped._gr. ~ Treatment, scales = "free_x") +
  geom_col(show.legend = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "p-value", y = "Log-fold change", x = "Taxonomic family") +
  scale_fill_manual(values = c("#95D840", "#33638D", "#481567"), 
                    drop = FALSE,
                    labels=c('<= 0.017', '<= 0.0033', '<= 0.00033')) +
  geom_hline(yintercept = 0, colour = "firebrick4") +
  theme(plot.margin = unit(c(1,1,1,2), "cm"))
```

Figure 5 in the article :

``` r
ggarrange(p_16S, p_ITS, 
          common.legend = TRUE, legend = "bottom", 
          labels = c("(a)", "(b)"), ncol = 1)
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/03_fig5-1.png)<!-- -->

``` r
ggsave(file.path(article_fig_path, "fig5.pdf"), 
       width = 12, height = 15, units = "in")
```

# 7 Session info

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=French_Canada.utf8  LC_CTYPE=French_Canada.utf8   
    ## [3] LC_MONETARY=French_Canada.utf8 LC_NUMERIC=C                  
    ## [5] LC_TIME=French_Canada.utf8    
    ## 
    ## time zone: America/Toronto
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices datasets  utils     methods   base     
    ## 
    ## other attached packages:
    ##  [1] conflicted_1.2.0  lubridate_1.9.5   forcats_1.0.1     stringr_1.6.0    
    ##  [5] dplyr_1.2.1       purrr_1.2.2       readr_2.2.0       tidyr_1.3.2      
    ##  [9] tibble_3.3.1      tidyverse_2.0.0   doRNG_1.8.6.3     rngtools_1.5.2   
    ## [13] foreach_1.5.2     ANCOMBC_2.8.1     ggpubr_1.0.0      ggplot2_4.0.3    
    ## [17] wesanderson_0.3.7 phyloseq_1.50.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      rstudioapi_0.19.0       jsonlite_2.0.0         
    ##   [4] magrittr_2.0.5          TH.data_1.1-5           nloptr_2.2.1           
    ##   [7] farver_2.1.2            rmarkdown_2.31          ragg_1.5.2             
    ##  [10] fs_2.1.0                zlibbioc_1.52.0         vctrs_0.7.3            
    ##  [13] multtest_2.62.0         memoise_2.0.1           minqa_1.2.8            
    ##  [16] base64enc_0.1-6         rstatix_1.1.0           htmltools_0.5.9        
    ##  [19] energy_1.7-12           haven_2.5.5             broom_1.0.13           
    ##  [22] cellranger_1.1.0        Rhdf5lib_1.28.0         Formula_1.2-6          
    ##  [25] rhdf5_2.50.2            htmlwidgets_1.6.4       sandwich_3.1-3         
    ##  [28] plyr_1.8.9              cachem_1.1.0            zoo_1.9-0              
    ##  [31] rootSolve_1.8.2.4       igraph_2.3.3            lifecycle_1.0.5        
    ##  [34] iterators_1.0.14        pkgconfig_2.0.3         Matrix_1.7-6           
    ##  [37] R6_2.6.1                fastmap_1.2.0           GenomeInfoDbData_1.2.13
    ##  [40] rbibutils_2.4.1         numDeriv_2016.8-1.1     digest_0.6.39          
    ##  [43] Exact_3.3               colorspace_2.1-3        S4Vectors_0.44.0       
    ##  [46] textshaping_1.0.5       Hmisc_5.2-6             vegan_2.7-5            
    ##  [49] labeling_0.4.3          timechange_0.4.0        httr_1.4.8             
    ##  [52] abind_1.4-8             mgcv_1.9-4              compiler_4.4.1         
    ##  [55] proxy_0.4-29            gsl_2.1-8               bit64_4.8.4            
    ##  [58] withr_3.0.3             doParallel_1.0.17       htmlTable_2.5.0        
    ##  [61] S7_0.2.2                backports_1.5.1         carData_3.0-6          
    ##  [64] ggsignif_0.6.4          MASS_7.3-66             biomformat_1.34.0      
    ##  [67] gtools_3.9.5            permute_0.9-10          CVXR_1.0-15            
    ##  [70] gld_2.6.8               tools_4.4.1             foreign_0.8-91         
    ##  [73] otel_0.2.0              ape_5.8-1               nnet_7.3-21            
    ##  [76] glue_1.8.1              nlme_3.1-170            rhdf5filters_1.18.1    
    ##  [79] grid_4.4.1              checkmate_2.3.4         cluster_2.1.8.3        
    ##  [82] reshape2_1.4.5          ade4_1.7-24             generics_0.1.4         
    ##  [85] gtable_0.3.6            tzdb_0.5.0              class_7.3-24           
    ##  [88] data.table_1.18.6.1     lmom_3.3                hms_1.1.4              
    ##  [91] car_3.1-5               XVector_0.46.0          BiocGenerics_0.52.0    
    ##  [94] pillar_1.11.1           splines_4.4.1           lattice_0.23-1         
    ##  [97] renv_1.1.4              survival_3.8-11         gmp_0.7-5.1            
    ## [100] bit_4.6.0               tidyselect_1.2.1        Biostrings_2.74.1      
    ## [103] knitr_1.51              reformulas_0.4.4        gridExtra_2.3.1        
    ## [106] IRanges_2.40.1          stats4_4.4.1            xfun_0.60              
    ## [109] expm_1.0-1              Biobase_2.66.0          skimr_2.2.2            
    ## [112] stringi_1.8.9           UCSC.utils_1.2.0        yaml_2.3.12            
    ## [115] boot_1.3-32             evaluate_1.0.5          codetools_0.2-20       
    ## [118] BiocManager_1.30.27     cli_3.6.6               rpart_4.1.27           
    ## [121] systemfonts_1.3.2       DescTools_0.99.60       Rdpack_2.6.6           
    ## [124] repr_1.1.7              Rcpp_1.1.2              GenomeInfoDb_1.42.3    
    ## [127] readxl_1.5.0            parallel_4.4.1          lme4_2.0-6             
    ## [130] Rmpfr_1.1-2             mvtnorm_1.4-2           lmerTest_3.2-1         
    ## [133] scales_1.4.0            e1071_1.7-17            crayon_1.5.3           
    ## [136] rlang_1.3.0             cowplot_1.2.0           multcomp_1.4-32
