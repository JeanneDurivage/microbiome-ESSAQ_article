PERMANOVA of climate and soil variables on pedologic groups
================

- [1 Starting data](#1-starting-data)
- [2 Previsualisation with a PCA](#2-previsualisation-with-a-pca)
- [3 PERMANOVA on pedologic groups](#3-permanova-on-pedologic-groups)
- [4 Investigation of the significant results with permutational
  ANOVAs](#4-investigation-of-the-significant-results-with-permutational-anovas)
  - [4.1 Exploration of significant permutational
    ANOVAs](#41-exploration-of-significant-permutational-anovas)
- [5 Session info](#5-session-info)

``` r
library(FactoMineR) # PCoA
library(factoextra) # PCoA
```

    ## Le chargement a nécessité le package : ggplot2

    ## Warning: le package 'ggplot2' a été compilé avec la version R 4.4.3

    ## Welcome to factoextra!

    ## Want to learn more? See two factoextra-related books at https://www.datanovia.com/library/principal-component-methods

``` r
library(multcompView)   # Multiple comparaison results ( graph letters)
library(broom) # for tidy
library(rstatix)
```

    ## 
    ## Attachement du package : 'rstatix'

    ## L'objet suivant est masqué depuis 'package:stats':
    ## 
    ##     filter

``` r
library(RVAideMemoire) # permutational ANOVAs
```

    ## Warning: le package 'RVAideMemoire' a été compilé avec la version R 4.4.3

    ## *** Package RVAideMemoire v 0.9-83-12 ***

    ## 
    ## Attachement du package : 'RVAideMemoire'

    ## L'objet suivant est masqué depuis 'package:broom':
    ## 
    ##     bootstrap

``` r
library(coin) # permutational chi2
```

    ## Le chargement a nécessité le package : survival

    ## 
    ## Attachement du package : 'coin'

    ## Les objets suivants sont masqués depuis 'package:rstatix':
    ## 
    ##     chisq_test, conover_test, fligner_test, friedman_test,
    ##     kruskal_test, sign_test, wilcox_test

``` r
library(vegan) # PERMANOVA
```

    ## Le chargement a nécessité le package : permute

    ## Warning: le package 'permute' a été compilé avec la version R 4.4.3

``` r
library(microeco) 
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
theme_set(theme_bw()) # default ggplot theme

library(conflicted) # conflicts gestion
conflicts_prefer(dplyr::filter) 
```

    ## [conflicted] Will prefer dplyr::filter over any other package.

``` r
conflicts_prefer(dplyr::select) 
```

    ## [conflicted] Will prefer dplyr::select over any other package.

``` r
set.seed(146) # random.org
```

``` r
github_path       <- "C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article" # Change me!
fig_path          <- file.path(github_path, "figures")
article_fig_path  <- file.path(github_path, "figures_article")
data_path         <- file.path(github_path, "data")
private_data_path <- file.path(github_path, "data_private")

knitr::opts_chunk$set(fig.path = paste0(fig_path, "/"))
```

# 1 Starting data

For legal reasons, the raw data is not publically available.

``` r
# Import datasets
otu_16S   <- read.csv(file.path(private_data_path, "otu_16S.csv"), row.names = 1)
otu_ITS   <- read.csv(file.path(private_data_path, "otu_ITS.csv"), row.names = 1)
otu_C_fct <- read.csv(file.path(private_data_path, "otu_C_functions.csv"), row.names = 1)
env_raw   <- read.csv(file.path(private_data_path, "environment_data.csv"), row.names = 1)
```

``` r
# dbRDA-associated analyses from this script
dbRDA_16S   <- read.csv(file.path(data_path, "dbRDA_16S_results.csv"))
dbRDA_ITS   <- read.csv(file.path(data_path, "dbRDA_ITS_results.csv"))
dbRDA_C_fct <- read.csv(file.path(data_path, "dbRDA_C_fct_results.csv"))
```

We will convert character vectors into ordered factors.

``` r
env_raw <- env_raw |>
  mutate(
    Ped._gr. = factor(Ped._gr., levels = c(
      "Humo-ferric podzol", "Orthic gleysol",
      "Humic orthic gleysol"
    )),
    OF_type = factor(OF_type, levels = c(
      "No-OF", "L-Hog", "L-Dairy",
      "Other-OF", "S-Poultry"
    )),
    Crop = factor(Crop, levels = c("Maïs", "Soya", "Céréale", "Prairie"))
    # Corn, Soybean, Cereals, Prairie
  )
```

Extract variables names significant in \>= 1 dbRDA (see the
`01_biogeographic_drivers_analysis` script).

``` r
signif_vars <-
  rbind(dbRDA_16S, dbRDA_C_fct, dbRDA_ITS) |>
  filter(p_value <= 0.05 &
    !variable %in% c(
      "CDI", "Crop", "Compaction_index",
      "Perennials", "STIR", "OF_type"
    )) |>
  pull(variable) |>
  unique()
```

The PERMANOVA analysis will be conducted on the following variables :

``` r
print(signif_vars)
```

    ## [1] "pH"           "Ped._gr."     "Clay_prc"     "Silt_prc"     "GDD"         
    ## [6] "Phosphorus"   "C_to_N_ratio"

Creating a dataset only with the desired variables :

``` r
env <- env_raw |> select(all_of(signif_vars)) |> 
                  drop_na()
```

Scaling response variables :

``` r
data <- env |> mutate(across(where(is.numeric), scale),
                      across(where(is.numeric), as.numeric))
```

# 2 Previsualisation with a PCA

``` r
data_num <- data |> mutate(across(everything(), as.numeric)) |> select(-Ped._gr.)
pca <- prcomp(data_num)
```

``` r
ggbiplot::ggbiplot(pca,
  groups = data$Ped._gr.,
  var.factor = 1,
  ellipse = TRUE, ellipse.level = 0.5, ellipse.alpha = 0.1,
  ellipse.linewidth = .5, ellipse.prob = 0.95,
  varname.size = 3,
  point.size = 2
) +
  labs(fill = "Pedologic group", color = "Pedologic group")
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/02_PCA_soil_climate-1.png)<!-- -->

# 3 PERMANOVA on pedologic groups

``` r
response <- data_num
```

``` r
perm <- adonis2(
  response ~ Ped._gr.,
  data = data,
  method = "euclidean",
  permutations = 9999
)
```

``` r
print(perm)
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = response ~ Ped._gr., data = data, permutations = 9999, method = "euclidean")
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2    94.23 0.19881 9.5533  1e-04 ***
    ## Residual 77   379.77 0.80119                  
    ## Total    79   474.00 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Application condition verification (heteroskedasticity)

``` r
distX     <- dist(response)
disp      <- betadisper(distX, data$Ped._gr.)
disp_test <- permutest(disp, permutations = 9999)

print(disp_test)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq     F N.Perm Pr(>F)
    ## Groups     2  0.224 0.11190 0.207   9999 0.8112
    ## Residuals 77 41.621 0.54053

# 4 Investigation of the significant results with permutational ANOVAs

Permutational ANOVAs on each response variable.

``` r
anova_res <- data.frame()

for (var in colnames(response)) {
  
  form <- stats::reformulate("Ped._gr.", var)
  mod <- perm.anova(formula = form, data = data, nperm = 9999, progress = FALSE)
  res <- rstatix::tidy(mod) |>
    dplyr::filter(term == "Ped._gr.") |>
    mutate(variable = var)

  anova_res <- anova_res |> rbind(res)
}
```

P-value adjustment for multiple testing.

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

knitr::kable(anova_res, digits = c(0, 0, 2, 0, 1, 1, 4, 4, 0))
```

| variable     | term      | sumsq |  df | meansq | statistic | p.value |  p_adj | signif_adj |
|:-------------|:----------|------:|----:|-------:|----------:|--------:|-------:|:-----------|
| GDD          | Ped.\_gr. | 18.51 |   2 |    9.3 |      11.8 |  0.0001 | 0.0006 | \*\*\*     |
| Clay_prc     | Ped.\_gr. | 18.84 |   2 |    9.4 |      12.1 |  0.0002 | 0.0010 | \*\*\*     |
| pH           | Ped.\_gr. | 14.67 |   2 |    7.3 |       8.8 |  0.0003 | 0.0012 | \*\*       |
| Silt_prc     | Ped.\_gr. | 16.35 |   2 |    8.2 |      10.0 |  0.0003 | 0.0012 | \*\*       |
| Phosphorus   | Ped.\_gr. | 14.57 |   2 |    7.3 |       8.7 |  0.0004 | 0.0012 | \*\*       |
| C_to_N_ratio | Ped.\_gr. | 11.29 |   2 |    5.6 |       6.4 |  0.0026 | 0.0026 | \*\*       |

Extracting significant response variables.

``` r
signif_vars <- anova_res |> filter(p_adj < 0.05) |> pull(variable)
```

## 4.1 Exploration of significant permutational ANOVAs

With permutational t-tests.

``` r
data_letters <- data.frame()

for (var in signif_vars) {
  posthoc <- pairwise.perm.t.test(
    data[[var]],
    data$Ped._gr.,
    nperm = 9999,
    p.method = "holm",
    progress = FALSE
  )

  p_tri <- posthoc$p.value
  gr <- levels(data$Ped._gr.)
  p_mat <- matrix(NA,
    nrow = length(gr), ncol = length(gr),
    dimnames = list(gr, gr)
  )
  p_mat[rownames(p_tri), colnames(p_tri)] <- p_tri
  p_mat[colnames(p_tri), rownames(p_tri)] <- t(p_tri)
  diag(p_mat) <- 1

  letters <- multcompView::multcompLetters(p_mat)$Letters

  data_letters_u <- data.frame(
    variable = var,
    Ped._gr. = names(letters),
    Letters = letters,
    row.names = NULL
  )


  data_letters <- rbind(data_letters, data_letters_u)
}
```

Data preparation to plotting results

``` r
# Long df
env_long <- env |>
  pivot_longer(cols = -Ped._gr., names_to = "variable") |>
  mutate(
    variable = gsub("_", " ", variable),
    variable = case_match(variable,
      "GDD" ~ "GDD, base 5°C (°C)",
      "Iron" ~ "Iron (mg / kg)",
      .default = gsub("prc", "%", variable)
    ),
    Ped._gr. = factor(Ped._gr., levels = c(
      "Humo-ferric podzol",
      "Orthic gleysol",
      "Humic orthic gleysol"
    ))
  )
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `variable = case_match(...)`.
    ## Caused by warning:
    ## ! `case_match()` was deprecated in dplyr 1.2.0.
    ## ℹ Please use `recode_values()` instead.

``` r
# add groups letters to data
letters_plot <- env_long |>
  group_by(variable, Ped._gr.) |>
  summarise(y_pos = max(value, na.rm = TRUE) * 1.1, .groups = "drop") |>
  left_join(
    data_letters |>
      mutate(
        variable = gsub("_", " ", variable),
        variable = case_match(variable,
          "GDD" ~ "GDD, base 5°C (°C)",
          "Iron" ~ "Iron (mg / kg)",
          .default = gsub("prc", "%", variable)
        )
      ),
    by = c("variable", "Ped._gr.")
  )


# permutational ANOVA results
permANOVA_res <- anova_res |>
  select(variable, p_adj) |>
  mutate(
    variable = gsub("_", " ", variable),
    variable = case_match(variable,
      "GDD" ~ "GDD, base 5°C (°C)",
      "Iron" ~ "Iron (mg / kg)",
      .default = gsub("prc", "%", variable)
    ),
    label = paste0("adj p = ", formatC(p_adj, digits = 2))
  )

# extra padding on the y axis
y_space <- env_long |>
  group_by(variable) |>
  summarise(
    y_min = 0,
    y_max = max(value, na.rm = TRUE) * 1.25
  )
```

Results plots. This corresponds to Fig. 4 in the article.

``` r
ggplot(env_long, aes(Ped._gr., value, fill = Ped._gr.)) +
  geom_boxplot(alpha = 0.8, width = 0.7) +
  geom_text(
    data = letters_plot, aes(x = Ped._gr., y = y_pos, label = Letters),
    fontface = "bold", size = 5
  ) +
  geom_text(
    data = permANOVA_res,
    aes(label = label),
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.2,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_blank(data = y_space, aes(y = y_max), inherit.aes = FALSE) +
  # geom_blank(data = y_space, aes(y = y_min), inherit.aes = FALSE) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_manual(values = c(
    "Humo-ferric podzol" = "#EF7924",
    "Orthic gleysol" = "#3480A2",
    "Humic orthic gleysol" = "#1C3C61"
  )) +
  labs(x = "Pedologic group", y = "Measured value", fill = " ") +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    plot.margin = margin(1, 1, 1, 1, unit = "cm")
  )
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/02_Fig4-1.png)<!-- -->

``` r
ggsave(file.path(article_fig_path, "fig4.pdf"), width = 15, height = 10, units = "in")
```

# 5 Session info

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
    ## [10] tidyverse_2.0.0         microeco_2.3.0          vegan_2.7-5            
    ## [13] permute_0.9-10          coin_1.4-5              survival_3.8-11        
    ## [16] RVAideMemoire_0.9-83-12 rstatix_1.1.0           broom_1.0.13           
    ## [19] multcompView_0.1-12     factoextra_2.2.0        ggplot2_4.0.3          
    ## [22] FactoMineR_2.16        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sandwich_3.1-3       rlang_1.3.0          magrittr_2.0.5      
    ##  [4] multcomp_1.4-32      otel_0.2.0           matrixStats_1.5.0   
    ##  [7] compiler_4.4.1       mgcv_1.9-4           systemfonts_1.3.2   
    ## [10] vctrs_0.7.3          reshape2_1.4.5       sysfonts_0.8.9      
    ## [13] pkgconfig_2.0.3      fastmap_1.2.0        backports_1.5.1     
    ## [16] labeling_0.4.3       rmarkdown_2.31       tzdb_0.5.0          
    ## [19] ragg_1.5.2           xfun_0.60            modeltools_0.2-24   
    ## [22] cachem_1.1.0         showtext_0.9-8       flashClust_1.1-4    
    ## [25] irlba_2.3.7          parallel_4.4.1       cluster_2.1.8.3     
    ## [28] ggbiplot_0.6.2       R6_2.6.1             stringi_1.8.9       
    ## [31] RColorBrewer_1.1-3   car_3.1-5            estimability_2.0.0  
    ## [34] Rcpp_1.1.2           knitr_1.51           zoo_1.9-0           
    ## [37] timechange_0.4.0     Matrix_1.7-6         splines_4.4.1       
    ## [40] igraph_2.3.3         tidyselect_1.2.1     rstudioapi_0.19.0   
    ## [43] abind_1.4-8          yaml_2.3.12          ggtext_0.1.2        
    ## [46] codetools_0.2-20     lattice_0.23-1       plyr_1.8.9          
    ## [49] withr_3.0.3          S7_0.2.2             evaluate_1.0.5      
    ## [52] xml2_1.6.0           pillar_1.11.1        BiocManager_1.30.27 
    ## [55] carData_3.0-6        renv_1.1.4           DT_0.34.0           
    ## [58] stats4_4.4.1         generics_0.1.4       hms_1.1.4           
    ## [61] scales_1.4.0         xtable_1.8-8         leaps_3.2           
    ## [64] glue_1.8.1           emmeans_2.0.4        scatterplot3d_0.3-45
    ## [67] tools_4.4.1          data.table_1.18.6.1  mvtnorm_1.4-2       
    ## [70] grid_4.4.1           ape_5.8-1            libcoin_1.0-13      
    ## [73] nlme_3.1-170         showtextdb_3.0       Formula_1.2-6       
    ## [76] cli_3.6.6            textshaping_1.0.5    gtable_0.3.6        
    ## [79] digest_0.6.39        ggrepel_0.9.8        TH.data_1.1-5       
    ## [82] htmlwidgets_1.6.4    farver_2.1.2         memoise_2.0.1       
    ## [85] htmltools_0.5.9      lifecycle_1.0.5      gridtext_0.1.6      
    ## [88] MASS_7.3-66
