Differential abundance analysis of prokaryotic carbon cycling potential
================

- [1 Setup](#1-setup)
- [2 Data preparation](#2-data-preparation)
- [3 PERMANOVA](#3-permanova)
- [4 Pairwise PERMANOVAs](#4-pairwise-permanovas)
- [5 Permutational ANOVAs](#5-permutational-anovas)
- [6 Session info](#6-session-info)

# 1 Setup

This R projects uses `renv`. A `.lock` file to restore adequate package
versions is available as part of the GitHub. See the [package
vignette](https://rstudio.github.io/renv/articles/renv.html#collaboration)
to know more.

``` r
library(microeco)
```

    ## Le chargement a nécessité le package : ggplot2

    ## Warning: le package 'ggplot2' a été compilé avec la version R 4.4.3

``` r
library(vegan)          # For the PERMANOVA dans the aitchison distance
```

    ## Le chargement a nécessité le package : permute

    ## Warning: le package 'permute' a été compilé avec la version R 4.4.3

``` r
library(broom)          # For tidy
library(RVAideMemoire)  # For perm.anova
```

    ## Warning: le package 'RVAideMemoire' a été compilé avec la version R 4.4.3

    ## *** Package RVAideMemoire v 0.9-83-12 ***

    ## 
    ## Attachement du package : 'RVAideMemoire'

    ## L'objet suivant est masqué depuis 'package:broom':
    ## 
    ##     bootstrap

``` r
library(rstatix)
```

    ## 
    ## Attachement du package : 'rstatix'

    ## L'objet suivant est masqué depuis 'package:stats':
    ## 
    ##     filter

``` r
library(multcompView)   # Post-hoc groups - generatin group letters
library(pairwiseAdonis) # Pairwise PERMANOVAS
```

    ## Le chargement a nécessité le package : cluster

``` r
library(skimr)          # Summarising dataframes
```

    ## Warning: le package 'skimr' a été compilé avec la version R 4.4.3

``` r
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
    ## ✖ RVAideMemoire::bootstrap() masks broom::bootstrap()
    ## ✖ dplyr::filter()            masks rstatix::filter(), stats::filter()
    ## ✖ dplyr::lag()               masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
theme_set(theme_bw())

library(conflicted)
conflicts_prefer(dplyr::filter, dplyr::select)
```

    ## [conflicted] Will prefer dplyr::filter over any other package.
    ## [conflicted] Will prefer dplyr::select over any other package.

``` r
set.seed(389) # random.org
```

``` r
github_path       <- "C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article" # Change-me!
fig_path          <- file.path(github_path, "figures")
article_fig_path  <- file.path(github_path, "figures_article")
data_path         <- file.path(github_path, "data")
private_data_path <- file.path(github_path, "data_private")

knitr::opts_chunk$set(fig.path = paste0(fig_path, "/"))
```

# 2 Data preparation

For confidientality reasons, the raw data is not publically available.

``` r
otu_C_fct <- read.csv(file.path(private_data_path, "otu_C_functions.csv"), row.names = 1)
env_raw   <- read.csv(file.path(private_data_path, "environment_data.csv"), row.names = 1)
```

This scripts uses 2 datasets :

- `env_raw`, which contains all necessary soil-, climate- and
  farming-related variables.

- `otu_C_fct` which contains the selected and grouped PiCRUST2
  procaryotic functions prediction

We will prepare the `env` table by selecting necessary variables and
rows and scaling.

``` r
env <- env_raw |> 
  rownames_to_column("site") |> 
  filter(site %in% colnames(otu_C_fct)) |> 
  column_to_rownames("site") |> 
  select(OF_type, Ped._gr., STIR) |> 
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

Here is an overview of the environmental data :

``` r
skim(env)
```

|                                                  |      |
|:-------------------------------------------------|:-----|
| Name                                             | env  |
| Number of rows                                   | 77   |
| Number of columns                                | 3    |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |      |
| Column type frequency:                           |      |
| factor                                           | 2    |
| numeric                                          | 1    |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |      |
| Group variables                                  | None |

Data summary

**Variable type: factor**

| skim_variable | n_missing | complete_rate | ordered | n_unique | top_counts |
|:---|---:|---:|:---|---:|:---|
| OF_type | 0 | 1 | FALSE | 5 | Oth: 22, No-: 21, L-D: 17, L-H: 10 |
| Ped.\_gr. | 0 | 1 | FALSE | 3 | Hum: 45, Hum: 20, Ort: 12 |

**Variable type: numeric**

| skim_variable | n_missing | complete_rate | mean |  sd |    p0 |   p25 |  p50 |  p75 | p100 | hist  |
|:--------------|----------:|--------------:|-----:|----:|------:|------:|-----:|-----:|-----:|:------|
| STIR          |         0 |             1 |    0 |   1 | -1.52 | -0.82 | 0.01 | 0.71 | 4.17 | ▇▇▅▁▁ |

Here is an overview of the procaryotic carbon cycling data :

``` r
reponse_data <- otu_C_fct |> t()

skim(reponse_data)
```

|                                                  |              |
|:-------------------------------------------------|:-------------|
| Name                                             | reponse_data |
| Number of rows                                   | 77           |
| Number of columns                                | 8            |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |              |
| Column type frequency:                           |              |
| numeric                                          | 8            |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |              |
| Group variables                                  | None         |

Data summary

**Variable type: numeric**

| skim_variable | n_missing | complete_rate | mean | sd | p0 | p25 | p50 | p75 | p100 | hist |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| Methanotrophy | 0 | 1 | 151751.57 | 42077.11 | 82965.96 | 121647.92 | 136612.67 | 180997.33 | 259997.37 | ▅▇▅▂▁ |
| Methanogenesis | 0 | 1 | 90015.55 | 28130.84 | 43637.72 | 65060.98 | 87014.93 | 111848.61 | 171541.19 | ▇▇▆▃▁ |
| Carbon_fixation | 0 | 1 | 462774.23 | 128509.85 | 256226.66 | 367231.45 | 421957.78 | 551957.33 | 795215.40 | ▅▇▃▃▁ |
| Labile_carbon_degradation | 0 | 1 | 411799.76 | 124492.23 | 211743.45 | 319873.87 | 381416.32 | 488853.89 | 745241.85 | ▅▇▃▂▂ |
| Recalcitrant_carbon_degradation | 0 | 1 | 612050.04 | 179513.65 | 325930.87 | 488149.41 | 568736.56 | 728621.00 | 1123971.51 | ▅▇▃▁▁ |
| Aerobic_respiration | 0 | 1 | 616919.00 | 172969.87 | 370329.72 | 503083.92 | 560328.58 | 734976.32 | 1048857.57 | ▆▇▂▂▂ |
| Anaerobic_respiration | 0 | 1 | 34816.51 | 10098.22 | 18042.73 | 27617.23 | 32076.44 | 41931.44 | 61679.47 | ▅▇▅▂▂ |
| Fermentation | 0 | 1 | 468243.23 | 127961.98 | 270072.14 | 378332.86 | 427501.87 | 561924.02 | 793275.68 | ▅▇▃▂▂ |

Joining all data in a single dataframe to simplify parts of the code :

``` r
data_all <- otu_C_fct |> rownames_to_column("C_function") |> 
                         pivot_longer(-C_function, names_to = "site") |> 
                         pivot_wider(names_from = C_function) |> 
                         left_join(env |> rownames_to_column("site")) |> 
                         column_to_rownames("site")
```

    ## Joining with `by = join_by(site)`

# 3 PERMANOVA

First, let’s see if there is an effect of STIR or OF type, individually
or in an interaction with the pedological group. The model will also
test wether the pedological group is significant by itself.

``` r
perm_global <- adonis2(reponse_data ~ Ped._gr.*(STIR + OF_type),
                       data = data_all,
                       permutations = 9999,
                       method = "aitchison")
```

``` r
print(perm_global)
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = reponse_data ~ Ped._gr. * (STIR + OF_type), data = data_all, permutations = 9999, method = "aitchison")
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model    16   1.8538 0.43618 2.9011  1e-04 ***
    ## Residual 60   2.3963 0.56382                  
    ## Total    76   4.2500 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Model is significant. Let’s see if it is explained by the interactions :

``` r
perm_interactions <- adonis2(reponse_data ~ Ped._gr.*(STIR + OF_type),
                             data = data_all,
                             permutations = 9999,
                             method = "aitchison",
                             by="margin")
```

``` r
print(perm_interactions)
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = reponse_data ~ Ped._gr. * (STIR + OF_type), data = data_all, permutations = 9999, method = "aitchison", by = "margin")
    ##                  Df SumOfSqs      R2      F Pr(>F)
    ## Ped._gr.:STIR     2   0.0250 0.00588 0.3129 0.9094
    ## Ped._gr.:OF_type  7   0.2354 0.05538 0.8419 0.6233
    ## Residual         60   2.3963 0.56382              
    ## Total            76   4.2500 1.00000

The interactions are not significant. Let’s now look at the individual
effects of OF type and STIR to determine if they are significant :

``` r
perm_single <- adonis2(reponse_data ~ Ped._gr.+ STIR + OF_type,
                       data = data_all,
                       permutations = 9999,
                       method = "aitchison",
                       by = "margin")
```

``` r
print(perm_single)
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = reponse_data ~ Ped._gr. + STIR + OF_type, data = data_all, permutations = 9999, method = "aitchison", by = "margin")
    ##          Df SumOfSqs      R2       F Pr(>F)    
    ## Ped._gr.  2   1.1548 0.27171 14.9762 0.0001 ***
    ## STIR      1   0.0238 0.00561  0.6184 0.5504    
    ## OF_type   4   0.3461 0.08142  2.2440 0.0277 *  
    ## Residual 69   2.6602 0.62593                   
    ## Total    76   4.2500 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

There is a significant effet of the OF type and the pedological group.

Verigication of the application conditions (we want it to be non
significant).

``` r
dist_X <- vegdist(reponse_data, "aitchison")

disp_pedo <- betadisper(dist_X, data_all$Ped._gr.)
print(anova(disp_pedo))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq F value  Pr(>F)  
    ## Groups     2 0.03956 0.0197784   2.449 0.09336 .
    ## Residuals 74 0.59764 0.0080762                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
disp_OF <- betadisper(dist_X, data_all$OF_type)
print(anova(disp_OF))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     4 0.03873 0.0096818  1.2101 0.3141
    ## Residuals 72 0.57608 0.0080012

The homogeneity of dispersion condition was respected.

# 4 Pairwise PERMANOVAs

1)  Pairwise PERMANOVAs on the **OF type** to determine which group(s)
    significantly differ from the other(s).

``` r
pairwise_OF <- pairwise.adonis2(
  reponse_data ~ OF_type,
  data = data_all,
  permutations = 9999,
  method = "aitchison"
)
```

Extracting p-values and R2 + p-value correction for multiple comparisons  

``` r
p_vals_OF <- sapply(pairwise_OF[-1], function(x) x$`Pr(>F)`[1])
R2_OF <- sapply(pairwise_OF[-1], function(x) x$R2[1])

data.frame(R2 = R2_OF, 
           p = p_vals_OF) |> 
  rownames_to_column("comparison") |> 
  filter(!grepl("Other", comparison)) |> 
  mutate(p.adj = p.adjust(p, method = "holm")) |> 
  print()
```

    ##             comparison         R2      p  p.adj
    ## 1   S-Poultry_vs_L-Hog 0.09885331 0.1885 0.7156
    ## 2   S-Poultry_vs_No-OF 0.04342413 0.2882 0.7156
    ## 3 S-Poultry_vs_L-Dairy 0.03829885 0.4282 0.7156
    ## 4       L-Hog_vs_No-OF 0.16914807 0.0075 0.0450
    ## 5     L-Hog_vs_L-Dairy 0.10079874 0.0607 0.3035
    ## 6     No-OF_vs_L-Dairy 0.04323282 0.1789 0.7156

There is only a marginal difference between `No-OF` and `L-Hog`.

2)  Pairwise PERMANOVAs on the **pedological group** to determine which
    group(s) significantly differ from the other(s).

``` r
pairwise_pedo <- pairwise.adonis2(reponse_data ~ Ped._gr.,
                                  data = data_all,
                                  permutations = 9999,
                                  method = "aitchison")
```

Extracting p-values and R2 + p-value correction for multiple comparisons  

``` r
p_vals_pedo <- sapply(pairwise_pedo[-1], function(x) x$`Pr(>F)`[1])
R2_pedo <- sapply(pairwise_pedo[-1], function(x) x$R2[1])

data.frame(R2 = R2_pedo, 
           p = p_vals_pedo) |> 
  rownames_to_column("comparison") |> 
  filter(!grepl("Other", comparison)) |> 
  mutate(p.adj = p.adjust(p, method = "holm")) |> 
  print()
```

    ##                                   comparison         R2      p  p.adj
    ## 1 Humo-ferric podzol_vs_Humic orthic gleysol 0.30735925 0.0001 0.0003
    ## 2       Humo-ferric podzol_vs_Orthic gleysol 0.23433543 0.0001 0.0003
    ## 3     Humic orthic gleysol_vs_Orthic gleysol 0.05289692 0.0424 0.0424

Humo-ferric podzols differ from both gleysols. There is a small, barely
significant difference between the two gleysol groups.

# 5 Permutational ANOVAs

Considering the important differences between pedological groups, let’s
investigate which carbon cycling functions differ between pedological
groups using permutational ANOVAs

Running permutational ANOVAs on each of the carbon cycling functions :

``` r
anova_res <- data.frame()

for(var in colnames(reponse_data)){

  var_nom <- var
  form <- reformulate("Ped._gr.", var_nom)
  mod <- perm.anova(formula = form, data = data_all, nperm = 9999, progress = FALSE)
  res <- tidy(mod) |> filter(term == "Ped._gr.") |> mutate(variable = var_nom)
  
  anova_res <- anova_res |> rbind(res)
}
```

P-value correction for multiple testing :

``` r
anova_res <- anova_res |> 
  mutate(p_adj = p.adjust(p.value, method = "holm"),
         signif_adj = case_when(
           p_adj <= 0.001 ~ "***",
           p_adj <= 0.01 ~ "**",
           p_adj <= 0.05 ~ "*",
           p_adj > 0.05 ~ "n.s."
         )) |> 
  arrange(p_adj) |> 
  relocate(variable, .before = term)
```

``` r
print(anova_res)
```

    ## # A tibble: 8 × 9
    ##   variable       term    sumsq    df  meansq statistic p.value  p_adj signif_adj
    ##   <chr>          <chr>   <dbl> <int>   <dbl>     <dbl>   <dbl>  <dbl> <chr>     
    ## 1 Labile_carbon… Ped.… 2.28e11     2 1.14e11      8.87  0.0002 0.0016 **        
    ## 2 Aerobic_respi… Ped.… 5.02e11     2 2.51e11     10.5   0.0003 0.0021 **        
    ## 3 Fermentation   Ped.… 2.36e11     2 1.18e11      8.66  0.0003 0.0021 **        
    ## 4 Methanotrophy  Ped.… 2.59e10     2 1.30e10      8.82  0.0007 0.003  **        
    ## 5 Recalcitrant_… Ped.… 4.87e11     2 2.44e11      9.19  0.0006 0.003  **        
    ## 6 Methanogenesis Ped.… 9.67e 9     2 4.84e 9      7.09  0.0018 0.0054 **        
    ## 7 Carbon_fixati… Ped.… 1.98e11     2 9.92e10      6.95  0.002  0.0054 **        
    ## 8 Anaerobic_res… Ped.… 8.75e 8     2 4.37e 8      4.71  0.0124 0.0124 *

All carbon cycling functions show significant differences between
pedological groups.

Now, for each function, let’s perform permutational pairwise comparisons
to identify between which pedological groups these functions differ.

``` r
data_letters <- data.frame()

for(var in colnames(reponse_data)){

  # post hoc permutational t-tests
  posthoc <- pairwise.perm.t.test(
    data_all[[var]],
    data_all$Ped._gr.,
    nperm = 9999,
    p.method = "holm",
    progress = FALSE
  )
  
  # creating a p-values matrix
  p_tri <- posthoc$p.value
  gr <- levels(data_all$Ped._gr.)
  p_mat <- matrix(NA, nrow = length(gr), ncol = length(gr),
                  dimnames = list(gr, gr))
  p_mat[rownames(p_tri), colnames(p_tri)] <- p_tri
  p_mat[colnames(p_tri), rownames(p_tri)] <- t(p_tri)
  diag(p_mat) <- 1
  
  # identifying distinct groups (assigning group letters)
  letters <- multcompView::multcompLetters(p_mat)$Letters
  
  
  # converting results in a dataframe for easier manipulation
  data_letters_u <- data.frame(variable = var,
                               Pedo = names(letters),
                               Letters = letters,
                               row.names = NULL)
  
  
  data_letters <- rbind(data_letters, data_letters_u)
}
```

Now, let’s plot these results.

The data needed to create these graph cannot be shared since it is part
of the raw data covered by a legal privacy agreement (total reads per
site, for every carbon cycling procaryotic function, created with the
sequencing results).

``` r
# Core plot data
data_long <- data_all |>
  select(-STIR, -OF_type) |>
  pivot_longer(cols = -c("Ped._gr."), names_to = "variable") |>
  mutate(variable_clean = gsub("_", " ", variable))

# Data for adding the letters (representing the post-hoc test groups) on top of the boxplots
letters_plot <- data_long |>
  group_by(variable, Ped._gr.) |>
  summarise(y_pos = max(value, na.rm = TRUE), .groups = "drop") |>
  left_join(data_letters |> rename(Ped._gr. = Pedo), by = c("variable", "Ped._gr.")) |>
  mutate(
    y_pos = y_pos * 1.1,
    Letters = if_else(is.na(Letters), "a", Letters),
    variable_clean = gsub("_", " ", variable)
  )

# Data for text annotations indicating the permutational anovas p-values
graph_pvals <- anova_res |> 
  select(variable, p_adj) |> 
  mutate(
    variable_clean = gsub("_", " ", variable),
    label = paste0("adj p = ", formatC(p_adj, digits = 2))
    )

# Data to add extra padding around the y axis
y_space <- data_long |>
  group_by(variable_clean) |>
  summarise(y_min = 0,
            y_max = max(value, na.rm = TRUE) * 1.25)
```

This is Figure 6 in the article :

``` r
data_long |> 
  ggplot(aes(Ped._gr., value, fill = Ped._gr.)) +
  geom_boxplot(alpha = 0.8, width = 0.7) +
  geom_text(
    data = letters_plot,
    aes(x = Ped._gr., y = y_pos, label = Letters),
    fontface = "bold",
    size = 4
  ) +
  geom_text(
    data = graph_pvals,
    aes(label = label),
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.2,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_blank(data = y_space, aes(y = y_max), inherit.aes = FALSE) +

  facet_wrap(~variable_clean, scales = "free_y") +
  
  scale_fill_manual(values = c(
    "Humo-ferric podzol"   = "#EF7924", 
    "Orthic gleysol"       = "#2480A2",               
    "Humic orthic gleysol" = "#2D3C61"
  )) +
  labs(
    x = "\nPedological group",
    y = "Number of reads",
    fill = "Pedological group"
  ) +

  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12),
    plot.margin = margin(1, 1, 1, 1, unit = "cm")
  )
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/04_fig6-1.png)<!-- -->

``` r
ggsave(file.path(article_fig_path, "fig6.pdf"), width = 15, height = 10, units = "in")
```

# 6 Session info

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
    ##  [1] conflicted_1.2.0        lubridate_1.9.5         forcats_1.0.1          
    ##  [4] stringr_1.6.0           dplyr_1.2.1             purrr_1.2.2            
    ##  [7] readr_2.2.0             tidyr_1.3.2             tibble_3.3.1           
    ## [10] tidyverse_2.0.0         skimr_2.2.2             pairwiseAdonis_0.4.1   
    ## [13] cluster_2.1.8.3         multcompView_0.1-12     rstatix_1.1.0          
    ## [16] RVAideMemoire_0.9-83-12 broom_1.0.13            vegan_2.7-5            
    ## [19] permute_0.9-10          microeco_2.3.0          ggplot2_4.0.3          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        xfun_0.60           lattice_0.23-1     
    ##  [4] tzdb_0.5.0          vctrs_0.7.3         tools_4.4.1        
    ##  [7] generics_0.1.4      parallel_4.4.1      pkgconfig_2.0.3    
    ## [10] Matrix_1.7-6        data.table_1.18.6.1 RColorBrewer_1.1-3 
    ## [13] S7_0.2.2            lifecycle_1.0.5     compiler_4.4.1     
    ## [16] farver_2.1.2        textshaping_1.0.5   repr_1.1.7         
    ## [19] carData_3.0-6       htmltools_0.5.9     yaml_2.3.12        
    ## [22] Formula_1.2-6       pillar_1.11.1       car_3.1-5          
    ## [25] MASS_7.3-66         cachem_1.1.0        abind_1.4-8        
    ## [28] nlme_3.1-170        tidyselect_1.2.1    digest_0.6.39      
    ## [31] stringi_1.8.9       reshape2_1.4.5      labeling_0.4.3     
    ## [34] splines_4.4.1       fastmap_1.2.0       grid_4.4.1         
    ## [37] cli_3.6.6           magrittr_2.0.5      base64enc_0.1-6    
    ## [40] utf8_1.2.6          ape_5.8-1           withr_3.0.3        
    ## [43] scales_1.4.0        backports_1.5.1     timechange_0.4.0   
    ## [46] rmarkdown_2.31      igraph_2.3.3        otel_0.2.0         
    ## [49] ragg_1.5.2          hms_1.1.4           memoise_2.0.1      
    ## [52] evaluate_1.0.5      knitr_1.51          mgcv_1.9-4         
    ## [55] rlang_1.3.0         Rcpp_1.1.2          glue_1.8.1         
    ## [58] BiocManager_1.30.27 renv_1.1.4          rstudioapi_0.19.0  
    ## [61] jsonlite_2.0.0      R6_2.6.1            plyr_1.8.9         
    ## [64] systemfonts_1.3.2
