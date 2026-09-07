dbRDAs - Relative importance of explanatory variables
================

- [1 Setup](#1-setup)
- [2 Starting data](#2-starting-data)
- [3 (Multi)collinearity diagnostic](#3-multicollinearity-diagnostic)
- [4 dbRDAs](#4-dbrdas)
  - [4.1 16S](#41-16s)
  - [4.2 ITS](#42-its)
  - [4.3 Carbon cycling functions](#43-carbon-cycling-functions)
  - [4.4 Final plot](#44-final-plot)
- [5 Variance partition analyses](#5-variance-partition-analyses)
  - [5.1 16S (ASV)](#51-16s-asv)
  - [5.2 ITS (ASV)](#52-its-asv)
  - [5.3 Cabon cycling functions](#53-cabon-cycling-functions)
  - [5.4 Final plot](#54-final-plot)
- [6 Session info](#6-session-info)

**This code shows the dbRDA analyses and the variance partitionning
analyses.**

# 1 Setup

``` r
set.seed(6) # random.org
```

This R projects uses `renv`. A `.lock` file to restore adequate package
versions is available as part of the GitHub. See the [package
vignette](https://rstudio.github.io/renv/articles/renv.html#collaboration)
to know more.

``` r
library(conflicted) # Package conflicts

library(knitr) # Table esthetics
```

    ## Warning: le package 'knitr' a été compilé avec la version R 4.4.3

``` r
library(ggpubr) # Puts multiple graphs together
```

    ## Le chargement a nécessité le package : ggplot2

    ## Warning: le package 'ggplot2' a été compilé avec la version R 4.4.3

``` r
library(patchwork) # Puts multiple graphs together
```

    ## Warning: le package 'patchwork' a été compilé avec la version R 4.4.3

``` r
library(wesanderson) # Graphs colours
```

    ## Warning: le package 'wesanderson' a été compilé avec la version R 4.4.3

``` r
library(microeco) # Sequencing data manipulation
library(vegan) # distances and statistical analyses
```

    ## Le chargement a nécessité le package : permute

    ## Warning: le package 'permute' a été compilé avec la version R 4.4.3

``` r
library(skimr) # Data summarisation
```

    ## Warning: le package 'skimr' a été compilé avec la version R 4.4.3

``` r
library(tidyverse) # Graphs + data manipulation
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

``` r
theme_set(theme_bw()) # Default ggplot theme
conflicts_prefer(dplyr::filter, dplyr::select) # Prevents conflicts between packages
```

    ## [conflicted] Will prefer dplyr::filter over any other package.

    ## [conflicted] Will prefer dplyr::select over any other package.

``` r
github_path       <- "C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article" # Change me!
fig_path          <- file.path(github_path, "figures/")
article_fig_path  <- file.path(github_path, "figures_article")
data_path         <- file.path(github_path, "data")
private_data_path <- file.path(github_path, "data_private")

knitr::opts_chunk$set(fig.path = paste0(fig_path, "/"))
```

# 2 Starting data

For confidientality reasons, the raw data is not publically available.

``` r
# Import datasets
otu_16S   <- read.csv(file.path(private_data_path, "otu_16S.csv"), row.names = 1) 
otu_ITS   <- read.csv(file.path(private_data_path, "otu_ITS.csv"), row.names = 1)
otu_C_fct <- read.csv(file.path(private_data_path, "otu_C_functions.csv"), row.names = 1)
env_raw   <- read.csv(file.path(private_data_path, "environment_data.csv"), row.names = 1) 
```

This scripts uses 4 datasets :

- `env`, which contains all soil-, climate- and farming-related
  variables.

- `otu_16S` and `otu_ITS`, which contains the ASV abundance table from
  the 16S and ITS rRNA amplicon sequencing

- `otu_C_fct`, which contains carbon cycling functionnal groups
  abundance.

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
    Crop = factor(Crop, levels = c("Maïs", "Soya", "Céréale", "Prairie")) # Corn, Soybean, Cereals, Prairie
  )
```

Here is an overview of the dataset :

``` r
skim(env_raw)
```

|                                                  |         |
|:-------------------------------------------------|:--------|
| Name                                             | env_raw |
| Number of rows                                   | 80      |
| Number of columns                                | 17      |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |         |
| Column type frequency:                           |         |
| factor                                           | 3       |
| numeric                                          | 14      |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |         |
| Group variables                                  | None    |

Data summary

**Variable type: factor**

| skim_variable | n_missing | complete_rate | ordered | n_unique | top_counts |
|:---|---:|---:|:---|---:|:---|
| Ped.\_gr. | 0 | 1 | FALSE | 3 | Hum: 46, Hum: 21, Ort: 13 |
| OF_type | 0 | 1 | FALSE | 5 | No-: 22, Oth: 22, L-D: 17, L-H: 12 |
| Crop | 0 | 1 | FALSE | 4 | Maï: 43, Soy: 22, Cér: 9, Pra: 6 |

**Variable type: numeric**

| skim_variable | n_missing | complete_rate | mean | sd | p0 | p25 | p50 | p75 | p100 | hist |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| pH | 0 | 1 | 6.60 | 0.55 | 5.48 | 6.26 | 6.58 | 6.92 | 7.84 | ▃▇▇▅▃ |
| Phosphorus | 0 | 1 | 75.49 | 68.52 | 9.19 | 29.54 | 55.09 | 101.33 | 418.74 | ▇▂▁▁▁ |
| Iron | 0 | 1 | 249.05 | 75.77 | 116.61 | 199.09 | 234.85 | 290.06 | 473.17 | ▅▇▆▂▁ |
| Carbon_prc | 0 | 1 | 1.91 | 0.86 | 0.64 | 1.33 | 1.71 | 2.37 | 4.26 | ▆▇▅▂▂ |
| C_to_N_ratio | 0 | 1 | 12.91 | 1.88 | 9.52 | 11.38 | 12.89 | 14.14 | 17.86 | ▆▇▇▅▁ |
| Silt_prc | 0 | 1 | 30.89 | 17.30 | 2.64 | 18.66 | 29.48 | 39.04 | 68.02 | ▆▇▇▃▃ |
| Clay_prc | 0 | 1 | 29.56 | 22.48 | 6.51 | 12.69 | 20.41 | 35.19 | 75.63 | ▇▃▁▁▃ |
| VESS | 0 | 1 | 2.31 | 0.67 | 1.00 | 1.90 | 2.42 | 2.75 | 3.58 | ▅▇▇▇▅ |
| PPT | 0 | 1 | 733.48 | 50.83 | 597.97 | 705.28 | 722.48 | 770.89 | 832.01 | ▁▂▇▂▃ |
| GDD | 0 | 1 | 2084.03 | 178.78 | 1635.08 | 1883.27 | 2191.47 | 2222.87 | 2273.64 | ▁▃▁▂▇ |
| STIR | 0 | 1 | 73.04 | 48.18 | 0.15 | 33.75 | 74.79 | 107.33 | 274.04 | ▇▇▅▁▁ |
| Comp.\_index | 0 | 1 | 20.69 | 14.36 | 0.00 | 11.46 | 17.59 | 25.56 | 75.62 | ▇▆▂▁▁ |
| Perennials | 0 | 1 | 0.48 | 1.08 | 0.00 | 0.00 | 0.00 | 0.00 | 4.00 | ▇▁▁▁▁ |
| CDI | 0 | 1 | 3.61 | 2.29 | 1.00 | 2.00 | 3.00 | 5.00 | 9.00 | ▇▅▁▂▂ |

``` r
print( otu_16S[1:6, 1:3] )
```

    ##           site_1 site_2 site_3
    ## ASV_16S_1    433    467    145
    ## ASV_16S_2     50     39      0
    ## ASV_16S_3     24    432      0
    ## ASV_16S_4    780   1118    711
    ## ASV_16S_5    986    811   1032
    ## ASV_16S_6      0    216     56

``` r
print( otu_ITS[1:6, 1:3] )
```

    ##           site_1 site_2 site_3
    ## ASV_16S_1      0     16      0
    ## ASV_16S_2     72     61      2
    ## ASV_16S_3    822     14      0
    ## ASV_16S_4      0     18      0
    ## ASV_16S_5      0      0      0
    ## ASV_16S_6      0     57      0

``` r
print( otu_C_fct[1:6, 1:3] )
```

    ##                                    site_1   site_2    site_3
    ## Methanotrophy                    219757.6 215350.4 129563.16
    ## Methanogenesis                   109649.7 111848.6  64973.42
    ## Carbon_fixation                  672758.4 668116.6 397622.41
    ## Labile_carbon_degradation        652253.9 665633.6 364748.92
    ## Recalcitrant_carbon_degradation 1009982.7 946627.5 555638.05
    ## Aerobic_respiration              986357.3 949603.6 577050.53

Let’s scale all numeric explanatory variables :

``` r
env <- env_raw |> mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))
```

# 3 (Multi)collinearity diagnostic

To acess multicollinearity, we will evaluate :

- The correlation. It should be \< 0.9

- The Variance inflation factor (VIF). It should be \< 10.

- The condition index. It should be \< 30.

**Correlation**

``` r
cor_env <- env_raw |>
  mutate(
    across(where(is.character), as.factor),
    across(where(is.factor), as.numeric)
  ) |>
  cor()

corrplot::corrplot(cor_env,
  type = "upper", order = "hclust",
  tl.col = "black", tl.srt = 45
)
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_explanatory_var_correlations-1.png)<!-- -->

``` r
cor_env |>
  as.data.frame() |>
  rownames_to_column("var1") |>
  pivot_longer(-var1, values_to = "cor", names_to = "var2") |>
  filter(var1 > var2) |>
  arrange(desc(cor)) |>
  head() |>
  knitr::kable(
    digits = 3,
    col.names = c("Variable 1", "Variable 2", "Pearson's coefficient"),
    caption = "Top correlartion coefficients"
  )
```

| Variable 1 | Variable 2 | Pearson’s coefficient |
|:-----------|:-----------|----------------------:|
| Perennials | CDI        |                 0.748 |
| Ped.\_gr.  | GDD        |                 0.484 |
| Ped.\_gr.  | Clay_prc   |                 0.470 |
| pH         | Ped.\_gr.  |                 0.429 |
| VESS       | Iron       |                 0.395 |
| Crop       | CDI        |                 0.373 |

Top correlartion coefficients

``` r
# Run diagnostics
mod_vif <- lm(rep(1, nrow(env_raw)) ~ ., data = env_raw)
diagnostic <- olsrr::ols_coll_diag(mod_vif)
```

**VIF:**

``` r
arrange(diagnostic$vif_t, desc(VIF)) |>
  filter(VIF > 3) |>
  # select(-Tolerance)
  knitr::kable(digits = 2)
```

| Variables                     | Tolerance |  VIF |
|:------------------------------|----------:|-----:|
| GDD                           |      0.20 | 5.12 |
| Clay_prc                      |      0.21 | 4.71 |
| Ped.\_gr.Humic orthic gleysol |      0.24 | 4.15 |
| Carbon_prc                    |      0.25 | 4.08 |
| CDI                           |      0.25 | 4.04 |
| Perennials                    |      0.29 | 3.46 |
| Ped.\_gr.Orthic gleysol       |      0.31 | 3.21 |

**Condition index:**

``` r
diagnostic$eig_cindex |>
  arrange(desc(`Condition Index`)) |>
  mutate(across(everything(), ~ round(.x, 2))) |>
  head() |>
  knitr::kable()
```

| Eigenvalue | Condition Index | intercept | pH | Phosphorus | Iron | Carbon_prc | C_to_N_ratio | Silt_prc | Clay_prc | VESS | Ped.\_gr.Orthic gleysol | Ped.\_gr.Humic orthic gleysol | PPT | GDD | STIR | Comp.\_index | OF_typeL-Hog | OF_typeL-Dairy | OF_typeOther-OF | OF_typeS-Poultry | Perennials | CDI | CropSoya | CropCéréale | CropPrairie |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 225.29 | 0.99 | 0.06 | 0.00 | 0.03 | 0.03 | 0.04 | 0.02 | 0.02 | 0.09 | 0.01 | 0.01 | 0.76 | 0.62 | 0.01 | 0.03 | 0.02 | 0.07 | 0.02 | 0.01 | 0.01 | 0.05 | 0.00 | 0.01 | 0.00 |
| 0.00 | 76.28 | 0.00 | 0.31 | 0.00 | 0.01 | 0.30 | 0.21 | 0.04 | 0.13 | 0.00 | 0.00 | 0.00 | 0.05 | 0.37 | 0.08 | 0.01 | 0.00 | 0.03 | 0.01 | 0.00 | 0.02 | 0.00 | 0.02 | 0.01 | 0.04 |
| 0.00 | 60.21 | 0.00 | 0.60 | 0.01 | 0.04 | 0.09 | 0.04 | 0.15 | 0.15 | 0.02 | 0.00 | 0.03 | 0.14 | 0.01 | 0.01 | 0.00 | 0.08 | 0.05 | 0.01 | 0.01 | 0.00 | 0.01 | 0.00 | 0.01 | 0.01 |
| 0.01 | 44.20 | 0.00 | 0.02 | 0.07 | 0.12 | 0.08 | 0.67 | 0.01 | 0.24 | 0.02 | 0.07 | 0.01 | 0.05 | 0.00 | 0.00 | 0.01 | 0.01 | 0.02 | 0.01 | 0.02 | 0.05 | 0.10 | 0.02 | 0.00 | 0.02 |
| 0.02 | 23.71 | 0.00 | 0.00 | 0.26 | 0.61 | 0.00 | 0.02 | 0.05 | 0.09 | 0.46 | 0.02 | 0.01 | 0.00 | 0.00 | 0.01 | 0.08 | 0.01 | 0.00 | 0.03 | 0.01 | 0.01 | 0.02 | 0.01 | 0.01 | 0.00 |
| 0.04 | 19.41 | 0.00 | 0.00 | 0.02 | 0.14 | 0.16 | 0.01 | 0.07 | 0.13 | 0.18 | 0.41 | 0.52 | 0.00 | 0.00 | 0.01 | 0.01 | 0.00 | 0.00 | 0.00 | 0.01 | 0.04 | 0.19 | 0.00 | 0.03 | 0.08 |

The condition index is extremely high because of the intercept. Let’s
evaluate it again, with the scaled data, that we will also use in the
analyses.

``` r
mod_vif <- lm(rep(1, nrow(env)) ~ ., data = env)
diagnostic <- olsrr::ols_coll_diag(mod_vif)

diagnostic$eig_cindex |>
  arrange(desc(`Condition Index`)) |>
  mutate(across(everything(), ~ round(.x, 2))) |>
  head() |>
  knitr::kable()
```

| Eigenvalue | Condition Index | intercept | pH | Phosphorus | Iron | Carbon_prc | C_to_N_ratio | Silt_prc | Clay_prc | VESS | Ped.\_gr.Orthic gleysol | Ped.\_gr.Humic orthic gleysol | PPT | GDD | STIR | Comp.\_index | OF_typeL-Hog | OF_typeL-Dairy | OF_typeOther-OF | OF_typeS-Poultry | Perennials | CDI | CropSoya | CropCéréale | CropPrairie |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.04 | 10.10 | 0.92 | 0.01 | 0.06 | 0.12 | 0.11 | 0.00 | 0.04 | 0.10 | 0.05 | 0.50 | 0.71 | 0.01 | 0.01 | 0.00 | 0.01 | 0.04 | 0.08 | 0.08 | 0.04 | 0.01 | 0.09 | 0.00 | 0.05 | 0.11 |
| 0.09 | 6.74 | 0.00 | 0.04 | 0.26 | 0.24 | 0.42 | 0.42 | 0.12 | 0.58 | 0.09 | 0.11 | 0.07 | 0.01 | 0.12 | 0.02 | 0.04 | 0.12 | 0.06 | 0.00 | 0.14 | 0.02 | 0.22 | 0.00 | 0.09 | 0.03 |
| 0.10 | 6.32 | 0.04 | 0.01 | 0.00 | 0.12 | 0.00 | 0.05 | 0.15 | 0.11 | 0.15 | 0.05 | 0.12 | 0.24 | 0.33 | 0.01 | 0.01 | 0.25 | 0.42 | 0.34 | 0.03 | 0.09 | 0.07 | 0.03 | 0.03 | 0.00 |
| 0.13 | 5.54 | 0.02 | 0.05 | 0.03 | 0.05 | 0.23 | 0.03 | 0.05 | 0.03 | 0.04 | 0.03 | 0.04 | 0.05 | 0.32 | 0.02 | 0.01 | 0.08 | 0.03 | 0.14 | 0.04 | 0.30 | 0.25 | 0.14 | 0.09 | 0.01 |
| 0.20 | 4.54 | 0.00 | 0.02 | 0.00 | 0.02 | 0.00 | 0.00 | 0.02 | 0.04 | 0.04 | 0.00 | 0.00 | 0.29 | 0.14 | 0.01 | 0.28 | 0.05 | 0.06 | 0.06 | 0.00 | 0.37 | 0.14 | 0.02 | 0.05 | 0.03 |
| 0.28 | 3.81 | 0.00 | 0.05 | 0.18 | 0.10 | 0.03 | 0.11 | 0.14 | 0.00 | 0.13 | 0.05 | 0.00 | 0.04 | 0.02 | 0.08 | 0.02 | 0.01 | 0.02 | 0.06 | 0.14 | 0.00 | 0.05 | 0.14 | 0.02 | 0.03 |

# 4 dbRDAs

Every dbRDA will be run three times : for 16S, ITS, then carbon cycling
datasets.

To limit code redundancy, we will create two functions. `dbRDA_analysis`
will perform the data analysis and `dbRDA_plotting` will create the
analysis result plots that are included in the research article.

``` r
dbRDA_analysis <- function(otu, env) {
  # otu: OTU table with sample names in rows
  # env: Environmental table with sample names in rows. Data must be scaled

  # ----- dbRDA -------
  dist <- vegdist(t(otu), method = "robust.aitchison")
  dbRDA_res <- capscale(dist ~ ., data = env)

  # ---- Model interpretation ----
  # Adjusted R2 (in %)
  var_expl <- signif(RsquareAdj(dbRDA_res)$adj.r.squared*100, 2)

  # Global pvalue
  global_signif <- anova.cca(dbRDA_res, permutations = 9999)

  # Significant variables identification
  anova_vars <- anova.cca(dbRDA_res, by = "margin", permutations = 9999)
  # saving results
  signf_vars <- anova_vars |>
    select(`Pr(>F)`) |>
    as.data.frame() |>
    rownames_to_column("variable") |>
    rename(p_value = `Pr(>F)`)


  # ---- Effect sizes ------
  # individual dbRDA on every variable, one by one
  partial_r2 <-
    purrr::map_dbl(names(env), function(variable) {
      form <- as.formula(paste("dist ~", variable))
      mod <- capscale(form, data = env)
      RsquareAdj(mod)$adj.r.squared
    })


  # ----- Function output -------
  indiv_effects <-
    data.frame(
      variable = names(env),
      Partial_R2 = partial_r2 # ,
      # Partial_signf = partial_signf
    ) |>
    arrange(desc(Partial_R2)) |>
    left_join(signf_vars, by = join_by("variable")) |>
    mutate(
      signif = case_when(p_value < 0.001 ~ "< 0.001",
        p_value < 0.01 ~ "< 0.01",
        p_value < 0.05 ~ "< 0.05",
        .default = "n.s."
      ),
      signif = factor(signif, levels = c(
        "< 0.001", "< 0.01",
        "< 0.05", "n.s."
      )),
      variable = factor(variable, levels = colnames(env))
    )


  return(list(
    dbRDA = dbRDA_res,
    global_signif = global_signif,
    indiv_effects = indiv_effects,
    var_expl = var_expl
  ))
}
```

``` r
dbRDA_plotting <- function(indiv_effects, var_expl, title,
                           ASV = TRUE) {
  # indiv_effects: table from the dbRDA_analysis function
  # type: “Prokaryotes”, “Fungi” or "Carbon cycling" (for plotting)
  # ASV: Specify FALSE for carbon cycling functions

  # ----- Add info about data category -----
  clim_names <- c("PPT", "GDD")
  farm_names <- c(
    "Perennials", "CDI",
    "Crop", "OF type",
    "Comp. index", "STIR"
  )

  soil_names <- c(
    "Ped. gr.", "pH", "Silt %", "Clay %", "Iron",
    "Phosphorus", "C to N ratio", "Carbon %", "VESS"
  )

  y_order <- c(clim_names, farm_names, rev(soil_names))

  indiv_effects <- indiv_effects |>
    mutate(
      variable = gsub("_", " ", variable),
      variable = gsub("prc", "%", variable),
      variable = fct_relevel(variable, y_order),
      group = case_when(
        variable %in% clim_names ~ "Climate",
        variable %in% farm_names ~ "Farming practices",
        variable %in% soil_names ~ "Soil"
      )
    )


  # ----- This function generates the plot parts for each data category -----
  plot_groupe_fixed <- function(indiv_effects,
                                show_x_axis_label = FALSE,
                                ASV = TRUE,
                                legend_position = "none") {
    ggplot(indiv_effects, aes(x = Partial_R2, y = variable, fill = signif)) +
      geom_col(width = 0.8, show.legend = TRUE) +
      scale_x_continuous(limits = if (ASV) {
        c(-0.015, 0.4)
      } else {
        c(-0.015, 0.4)
      }) +
      geom_vline(xintercept = 0, color = "grey5", linewidth = .5) +
      scale_fill_manual(
        values = c(
          "< 0.05" = "#95D840",
          "< 0.01" = "#33638D",
          "< 0.001" = "#481567",
          "n.s." = "grey85"
        ),
        drop = FALSE
      ) +
      labs(
        x = ifelse(show_x_axis_label, "R2 (partial dbRDAs)", ""),
        fill = "p-value (ANOVA CCA)",
        y = "",
        title = unique(indiv_effects$group)
      ) +
      theme_minimal() +
      theme(legend.position = legend_position)
  }

  # Create 1 plot per variable category
  if (ASV == TRUE) {
    plots <- list(
      list(
        plot = plot_groupe_fixed(filter(indiv_effects, group == "Soil")),
        height = nrow(filter(indiv_effects, group == "Soil"))
      ),
      list(
        plot = plot_groupe_fixed(filter(indiv_effects, group == "Climate")),
        height = nrow(filter(indiv_effects, group == "Climate"))
      ),
      list(
        plot = plot_groupe_fixed(filter(indiv_effects, group == "Farming practices"),
          show_x_axis_label = TRUE
        ),
        height = nrow(filter(indiv_effects, group == "Farming practices"))
      )
    )
  } else {
    plots <- list(
      list(
        plot = plot_groupe_fixed(filter(indiv_effects, group == "Soil"), ASV = FALSE),
        height = nrow(filter(indiv_effects, group == "Soil"))
      ),
      list(
        plot = plot_groupe_fixed(filter(indiv_effects, group == "Climate"), ASV = FALSE),
        height = nrow(filter(indiv_effects, group == "Climate"))
      ),
      list(
        plot = plot_groupe_fixed(filter(indiv_effects, group == "Farming practices"),
          show_x_axis_label = TRUE, ASV = FALSE,
          legend_position = "bottom"
        ),
        height = nrow(filter(indiv_effects, group == "Farming practices"))
      )
    )
  }

  # ----- Plot assembly with patchwork -----
  var_expl_lab <- paste("Explains", var_expl, "% of the variance.")

  patch <-
    plots[[1]]$plot + plots[[2]]$plot + plots[[3]]$plot +
      plot_layout(
        heights = lapply(plots, `[[`, "height"),
        guides = "collect"
      ) +
      plot_annotation(
        title = title,
        subtitle = var_expl_lab,
        theme = theme(
          plot.title = element_text(size = 15),
          plot.caption = element_text(size = 13)
        )
      ) &
      theme_minimal() &
      theme(legend.position = "bottom")


  return(patch)
}
```

## 4.1 16S

``` r
dbRDA_16S <- dbRDA_analysis(otu_16S, env)
```

**Global model signifiance :**

``` r
print(dbRDA_16S$global_signif)
```

    ## Permutation test for capscale under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: capscale(formula = dist ~ pH + Phosphorus + Iron + Carbon_prc + C_to_N_ratio + Silt_prc + Clay_prc + VESS + Ped._gr. + PPT + GDD + STIR + Comp._index + OF_type + Perennials + CDI + Crop, data = env)
    ##          Df Variance      F Pr(>F)    
    ## Model    23   727.01 13.157  1e-04 ***
    ## Residual 56   134.53                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Adjusted R squared :

``` r
print(RsquareAdj(dbRDA_16S$dbRDA)$adj.r.squared)
```

    ## [1] 0.7797113

**dbRDA illusration**

``` r
plot(dbRDA_16S$dbRDA, display = c("sites", "bp"))
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_dbRDA_16S-1.png)<!-- -->

**Variables p-value and partial R2s**

``` r
dbRDA_16S$indiv_effects |> 
  mutate(Partial_R2 = round(Partial_R2, 3),
         variable = gsub("_", " ", variable),
         variable = gsub("prc", "%", variable)) |> 
  knitr::kable()
```

| variable     | Partial_R2 | p_value | signif   |
|:-------------|-----------:|--------:|:---------|
| pH           |      0.381 |  0.0001 | \< 0.001 |
| Ped. gr.     |      0.303 |  0.0001 | \< 0.001 |
| Clay %       |      0.259 |  0.0001 | \< 0.001 |
| GDD          |      0.148 |  0.2675 | n.s.     |
| PPT          |      0.137 |  0.1450 | n.s.     |
| C to N ratio |      0.120 |  0.1274 | n.s.     |
| Silt %       |      0.113 |  0.0106 | \< 0.05  |
| Perennials   |      0.097 |  0.3581 | n.s.     |
| Crop         |      0.078 |  0.3348 | n.s.     |
| Phosphorus   |      0.071 |  0.2380 | n.s.     |
| CDI          |      0.067 |  0.5856 | n.s.     |
| OF type      |      0.047 |  0.7005 | n.s.     |
| Iron         |      0.023 |  0.5903 | n.s.     |
| Comp. index  |      0.007 |  0.9626 | n.s.     |
| Carbon %     |      0.006 |  0.7390 | n.s.     |
| STIR         |      0.004 |  0.7470 | n.s.     |
| VESS         |     -0.004 |  0.5401 | n.s.     |

## 4.2 ITS

``` r
dbRDA_ITS <- dbRDA_analysis(otu_ITS, env)
```

**Global model signifiance :**

``` r
print(dbRDA_ITS$global_signif)
```

    ## Permutation test for capscale under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: capscale(formula = dist ~ pH + Phosphorus + Iron + Carbon_prc + C_to_N_ratio + Silt_prc + Clay_prc + VESS + Ped._gr. + PPT + GDD + STIR + Comp._index + OF_type + Perennials + CDI + Crop, data = env)
    ##          Df Variance      F Pr(>F)    
    ## Model    23   555.55 16.399  1e-04 ***
    ## Residual 56    82.48                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

**Adjusted R2**

``` r
print(RsquareAdj(dbRDA_ITS$dbRDA)$adj.r.squared)
```

    ## [1] 0.8176241

**Model illusration**

``` r
plot(dbRDA_ITS$dbRDA, display = c("sites", "bp"))
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_dbRDA_ITS-1.png)<!-- -->

**Variables p-value and partial R2s**

``` r
dbRDA_ITS$indiv_effects |> 
  mutate(Partial_R2 = round(Partial_R2, 3),
         variable = gsub("_", " ", variable), 
         variable = gsub("prc", "%", variable), 
         ) |> 
  knitr::kable()
```

| variable     | Partial_R2 | p_value | signif   |
|:-------------|-----------:|--------:|:---------|
| GDD          |      0.356 |  0.0001 | \< 0.001 |
| pH           |      0.292 |  0.0001 | \< 0.001 |
| Ped. gr.     |      0.259 |  0.0056 | \< 0.01  |
| Clay %       |      0.231 |  0.0059 | \< 0.01  |
| OF type      |      0.225 |  0.0049 | \< 0.01  |
| Carbon %     |      0.222 |  0.1872 | n.s.     |
| PPT          |      0.182 |  0.0952 | n.s.     |
| CDI          |      0.144 |  0.0200 | \< 0.05  |
| Phosphorus   |      0.137 |  0.0291 | \< 0.05  |
| Perennials   |      0.127 |  0.5347 | n.s.     |
| C to N ratio |      0.105 |  0.0330 | \< 0.05  |
| Silt %       |      0.063 |  0.1358 | n.s.     |
| Crop         |      0.040 |  0.4752 | n.s.     |
| STIR         |      0.039 |  0.1028 | n.s.     |
| Comp. index  |      0.026 |  0.4380 | n.s.     |
| Iron         |      0.001 |  0.3808 | n.s.     |
| VESS         |      0.000 |  0.2225 | n.s.     |

## 4.3 Carbon cycling functions

Because one site lacks functional data, it will have its own associated
`env` dataset.

``` r
env_fct <- env |> rownames_to_column("site") |> 
                  filter(site %in% colnames(otu_C_fct)) |> 
                  column_to_rownames("site")
```

``` r
dbRDA_C_fct <- dbRDA_analysis(otu_C_fct, env_fct)
```

**Global model signifiance :**

``` r
print(dbRDA_C_fct$global_signif)
```

    ## Permutation test for capscale under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: capscale(formula = dist ~ pH + Phosphorus + Iron + Carbon_prc + C_to_N_ratio + Silt_prc + Clay_prc + VESS + Ped._gr. + PPT + GDD + STIR + Comp._index + OF_type + Perennials + CDI + Crop, data = env)
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model    23   3.1514 6.6097  1e-04 ***
    ## Residual 53   1.0987                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

**Adjusted R2**

``` r
print(RsquareAdj(dbRDA_C_fct$dbRDA)$adj.r.squared)
```

    ## [1] 0.6293084

**Model illusration**

``` r
plot(dbRDA_C_fct$dbRDA, display = c("sites", "bp"))
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_dbRDA_C_fct-1.png)<!-- -->

**Variables p-value and partial R2s**

``` r
dbRDA_C_fct$indiv_effects |> 
  mutate(Partial_R2 = round(Partial_R2, 3),
         variable = gsub("_", " ", variable), 
         variable = gsub("prc", "%", variable), 
         ) |> 
  knitr::kable()
```

| variable     | Partial_R2 | p_value | signif   |
|:-------------|-----------:|--------:|:---------|
| pH           |      0.338 |  0.0001 | \< 0.001 |
| Clay %       |      0.299 |  0.0084 | \< 0.01  |
| Ped. gr.     |      0.267 |  0.0149 | \< 0.05  |
| PPT          |      0.142 |  0.2175 | n.s.     |
| Silt %       |      0.142 |  0.0216 | \< 0.05  |
| Phosphorus   |      0.136 |  0.8722 | n.s.     |
| C to N ratio |      0.128 |  0.4849 | n.s.     |
| GDD          |      0.095 |  0.7027 | n.s.     |
| OF type      |      0.040 |  0.2819 | n.s.     |
| Comp. index  |      0.038 |  0.1670 | n.s.     |
| Iron         |      0.016 |  0.7925 | n.s.     |
| VESS         |      0.010 |  0.1005 | n.s.     |
| Crop         |      0.008 |  0.2839 | n.s.     |
| Perennials   |     -0.001 |  0.0734 | n.s.     |
| STIR         |     -0.001 |  0.0725 | n.s.     |
| Carbon %     |     -0.004 |  0.4988 | n.s.     |
| CDI          |     -0.004 |  0.1116 | n.s.     |

## 4.4 Final plot

Saving the necessary data to plot results in the github :

``` r
# github_path <- "~/microbiome-ESSAQ_article" # Uncomment and change me
data_path   <- file.path(github_path, "data")

write.csv(dbRDA_16S$indiv_effects, file.path(data_path, "dbRDA_16S_results.csv"))
write.csv(dbRDA_ITS$indiv_effects, file.path(data_path, "dbRDA_ITS_results.csv"))
write.csv(dbRDA_C_fct$indiv_effects, file.path(data_path, "dbRDA_C_fct_results.csv"))
```

Load the data using this next chunk :

``` r
dbRDA_16S <- list()
dbRDA_ITS <- list()
dbRDA_C_fct <- list()

# dbRDA-associated analyses from this script
dbRDA_16S$indiv_effects   <- read.csv(file.path(data_path, "dbRDA_16S_results.csv"), 
                                      stringsAsFactors = TRUE)
dbRDA_ITS$indiv_effects   <- read.csv(file.path(data_path, "dbRDA_ITS_results.csv"), 
                                      stringsAsFactors = TRUE)
dbRDA_C_fct$indiv_effects <- read.csv(file.path(data_path, "dbRDA_C_fct_results.csv"), 
                                      stringsAsFactors = TRUE)

# % of variance explained by the global dbRDAs
dbRDA_16S$var_expl   <- 78
dbRDA_ITS$var_expl   <- 82
dbRDA_C_fct$var_expl <- 63
```

Creating individual plots using the pre-defined function :

``` r
# Creating the plots
patch_16S <- dbRDA_plotting(dbRDA_16S$indiv_effects,
                                 dbRDA_16S$var_expl,
                                 "(a) 16S ASVs")
patch_ITS <- dbRDA_plotting(dbRDA_ITS$indiv_effects,
                                 dbRDA_ITS$var_expl,
                                 "(c) ITS ASVs")
patch_C_fct <- dbRDA_plotting(dbRDA_C_fct$indiv_effects,
                                   dbRDA_C_fct$var_expl,
                                   "(b) Carbon cycling", 
                                   ASV = FALSE)
```

Assembling all plots to create Fig.2 from the article :

``` r
# tweaking some parameters to make it work...
patch_16S <- patch_16S & theme(plot.margin = margin(5, 5, 22, 5),
                               legend.position = "none",
                               axis.text.y = element_text(size = 12),
                               axis.text.x = element_text(size = 10)) 
patch_C_fct <- patch_C_fct & theme(
                               axis.text.y = element_text(size = 12),
                               axis.text.x = element_text(size = 10))

patch_ITS <- patch_ITS & theme(plot.margin = margin(5, 5, 22, 5),
                                   legend.position = "none",
                               axis.text.y = element_text(size = 12),
                               axis.text.x = element_text(size = 10))

# Assembling all plots
ggarrange(patch_16S, patch_C_fct, patch_ITS,
          nrow = 1)
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_Fig2-1.png)<!-- -->

``` r
ggsave("fig2.pdf", device = "pdf", path = article_fig_path, create.dir = TRUE,
       width = 13, height = 9, units = c("in"))
```

# 5 Variance partition analyses

We will use the exact same code organisation to conduct the variance
partition analyses.

``` r
varpart_analysis <- function(otu, env) {
  # ---- Data preparation (split data into categories) ----
  farm <- env |> select(
    "Perennials", "CDI", "Crop", "OF_type",
    "Comp._index", "STIR"
  )
  clim <- env |> select("PPT", "GDD")
  soil <- env |> select(
    "Phosphorus", "C_to_N_ratio", "Carbon_prc", "VESS",
    "Iron", "Silt_prc", "Clay_prc", "pH", "Ped._gr."
  )


  # ---- Variance partition analysis ----
  otu_rclr <- decostand(t(otu), method = "rclr")
  res_varpart <- varpart(otu_rclr, farm, soil, clim)

  # ---- Model interpretation ----
  otu_dist <- dist(otu_rclr)
  # Global significance
  dbRDA_total <- capscale(otu_dist ~ ., data = cbind(farm, soil, clim))
  signf_global <- anova.cca(dbRDA_total, permutations = 9999)

  # Significance of individual categories

  # farming variables
  vars_farm <- paste(names(farm), collapse = " + ")
  vars_soil_clim <- paste(c(names(soil), names(clim)), collapse = " + ")
  fml_farm <- as.formula(paste("otu_dist ~", vars_farm, "+ Condition(", vars_soil_clim, ")"))
  dbRDA_farm <- capscale(fml_farm, data = cbind(farm, soil, clim))
  signf_farm <- anova.cca(dbRDA_farm, permutations = 9999) # global test

  # soil variables
  vars_soil <- paste(names(soil), collapse = " + ")
  vars_farm_clim <- paste(c(names(farm), names(clim)), collapse = " + ")
  fml_soil <- as.formula(paste("otu_dist ~", vars_soil, "+ Condition(", vars_farm_clim, ")"))
  dbRDA_soil <- capscale(fml_soil, data = cbind(farm, soil, clim))
  signf_soil <- anova.cca(dbRDA_soil, permutations = 9999)

  # clim variables
  vars_clim <- paste(names(clim), collapse = " + ")
  vars_farm_soil <- paste(c(names(farm), names(soil)), collapse = " + ")
  fml_clim <- as.formula(paste("otu_dist ~", vars_clim, "+ Condition(", vars_farm_soil, ")"))
  dbRDA_clim <- capscale(fml_clim, data = cbind(farm, soil, clim))
  signf_clim <- anova.cca(dbRDA_clim, permutations = 9999)

  # Saving tests results in a list
  signifs <- data.frame(
    Group = c("Global model", "Farming practices", "Soil", "Climate"),
    model_p = c(
      signf_global$`Pr(>F)`[1],
      signf_farm$`Pr(>F)`[1],
      signf_soil$`Pr(>F)`[1],
      signf_clim$`Pr(>F)`[1]
    )
  ) |>
    mutate(model_padj = p.adjust(model_p))

  # ---- Total contributions per group ----
  res <- res_varpart$part$indfract |>
    as.data.frame() |>
    select(Adj.R.square) |>
    rownames_to_column("part")


  # Results formatting
  df <- tibble(
    group = c(
      "Farming practices", "Soil", "Climate",
      "Farming practices_X_Soil",
      "Farming practices_X_Climate",
      "Soil_X_Climate",
      "Farming practices_X_Soil_X_Climate", "Unexplained"
    ),
    value = res$Adj.R.square
  ) |>
    mutate(
      prc = value * 100,
      label = paste0(sprintf("%.1f", prc), "%"),
      label_legend = paste0(group, " (", sprintf("%.1f", prc), "%)")
    )


  df_groups <- list(
    `Farming practices` = df |> filter(str_detect(group, "Farming practices")),
    Soil                = df |> filter(str_detect(group, "Soil")),
    Climate             = df |> filter(str_detect(group, "Climate"))
  )

  groups_prop <-
    sapply(df_groups, function(x) sum(x$value)) |>
    as.data.frame() |>
    rename(Percentage = 1) |>
    mutate(Percentage = signif(Percentage * 100, digits = 3)) |>
    rownames_to_column(var = "Group")


  # ---- Output ----

  return(list(
    varpart = res_varpart,
    groups_prop = groups_prop,
    signifs = signifs,
    results_plot = df
  ))
}
```

``` r
plot_varpart <- function(results_plot, plot_title) {
  # ------ Legend settings ------
  groups_order <- c(
    "Soil",
    "Farming practices",
    "Climate",
    "Farming practices_X_Soil",
    "Soil_X_Climate",
    "Farming practices_X_Climate",
    "Farming practices_X_Soil_X_Climate",
    "Unexplained"
  )

  legend_colors <- c(
    "Soil"                               = "#DC3220",
    "Farming practices"                  = "#E1C901",
    "Climate"                            = "#1E88E5",
    "Farming practices_X_Soil"           = "#F28E2B",
    "Soil_X_Climate"                     = "#DB9A99",
    "Farming practices_X_Climate"        = "#A6CEE3",
    "Farming practices_X_Soil_X_Climate" = "purple4",
    "Unexplained"                        = "#DDDDDD"
  )


  group_labels <- list(
    Soil = "Soil",
    "Farming practices" = "Farming practices",
    Climate = "Climate",
    "Farming practices_X_Soil" = expression("Farming practices" %cap% "Soil"),
    "Soil_X_Climate" = expression("Soil" %cap% "Climate"),
    "Farming practices_X_Climate" = expression("Farming practices" %cap% "Climate"),
    "Farming practices_X_Soil_X_Climate" = expression("Farming practices" %cap% "Soil" %cap% "Climate"),
    Unexplained = "Unexplained"
  )

  group_labels <- list(
    Soil = "Soil",
    "Farming practices" = "Farming practices",
    Climate = "Climate",
    "Farming practices_X_Soil" = expression("Farming practices" ~ intersect() ~ "Soil"),
    "Soil_X_Climate" = expression("Soil" ~ intersect() ~ "Climate"),
    "Farming practices_X_Climate" =
      expression("Farming practices" ~ intersect() ~ "Climate"),
    "Farming practices_X_Soil_X_Climate" =
      expression("Farming practices" ~ intersect() ~ "Soil" ~ intersect() ~ "Climate"),
    Unexplained = "Unexplained"
  )


  # ------ Extract varpart results ------
  # res <- res_varpart$part$indfract$Adj.R.square

  # Results formatting
  df <- results_plot |> mutate(group = factor(group, levels = groups_order))

  # Donut plot
  p_part <-
    ggplot(df, aes(x = 2, y = value, fill = group)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    theme_void() +
    scale_fill_manual(
      values = legend_colors,
      breaks = names(legend_colors),
      labels = group_labels
    ) +
    labs(title = plot_title, fill = expression("Adjusted " * R^2))

  return(p_part)
}
```

For reference, here is how the variance part are named :

``` r
showvarparts(3)
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_reference_venn_plot-1.png)<!-- -->

## 5.1 16S (ASV)

**Variance partitionning :**

``` r
varpart_16S <- varpart_analysis(otu_16S, env)
```

**Variables groups + global significativity :**

| Group             | model_p | model_padj |
|:------------------|--------:|-----------:|
| Global model      |  0.0001 |     0.0004 |
| Farming practices |  0.2771 |     0.2771 |
| Soil              |  0.0001 |     0.0004 |
| Climate           |  0.0956 |     0.1912 |

**varpart() output**

``` r
print(varpart_16S$varpart)
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = otu_rclr, X = farm, soil, clim)
    ## 
    ## Explanatory tables:
    ## X1:  farm
    ## X2:  soil
    ## X3:  clim 
    ## 
    ## No. of explanatory tables: 3 
    ## Total variation (SS): 68062 
    ##             Variance: 861.54 
    ## No. of observations: 80 
    ## 
    ## Partition table:
    ##                       Df R.square Adj.R.square Testable
    ## [a+d+f+g] = X1        11  0.23529      0.11159     TRUE
    ## [b+d+e+g] = X2        10  0.78660      0.75567     TRUE
    ## [c+e+f+g] = X3         2  0.20866      0.18811     TRUE
    ## [a+b+d+e+f+g] = X1+X2 21  0.83323      0.77285     TRUE
    ## [a+c+d+e+f+g] = X1+X3 13  0.35027      0.22229     TRUE
    ## [b+c+d+e+f+g] = X2+X3 12  0.80798      0.77358     TRUE
    ## [a+b+c+d+e+f+g] = All 23  0.84385      0.77971     TRUE
    ## Individual fractions                                   
    ## [a] = X1 | X2+X3      11               0.00613     TRUE
    ## [b] = X2 | X1+X3      10               0.55742     TRUE
    ## [c] = X3 | X1+X2       2               0.00686     TRUE
    ## [d]                    0               0.02805    FALSE
    ## [e]                    0               0.10383    FALSE
    ## [f]                    0               0.01105    FALSE
    ## [g]                    0               0.06636    FALSE
    ## [h] = Residuals                        0.22029    FALSE
    ## Controlling 1 table X                                  
    ## [a+d] = X1 | X3       11               0.03418     TRUE
    ## [a+f] = X1 | X2       11               0.01718     TRUE
    ## [b+d] = X2 | X3       10               0.58547     TRUE
    ## [b+e] = X2 | X1       10               0.66126     TRUE
    ## [c+e] = X3 | X1        2               0.11069     TRUE
    ## [c+f] = X3 | X2        2               0.01791     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(varpart_16S$varpart, cutoff = 0, digits = 1, bg = c("#E1C901", "#DC3220", "#1E88E5"))
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_venn_plot_16S-1.png)<!-- -->

**TOTAL contribution per variable group**

| Group             | Percentage |
|:------------------|-----------:|
| Farming practices |       20.4 |
| Soil              |       66.3 |
| Climate           |       18.8 |

## 5.2 ITS (ASV)

**Variance partitionning :**

``` r
varpart_ITS <- varpart_analysis(otu_ITS, env)
```

**Variables groups + global significativity :**

``` r
knitr::kable(varpart_ITS$signifs)
```

| Group             | model_p | model_padj |
|:------------------|--------:|-----------:|
| Global model      |  0.0001 |     0.0004 |
| Farming practices |  0.0017 |     0.0017 |
| Soil              |  0.0001 |     0.0004 |
| Climate           |  0.0001 |     0.0004 |

**varpart() output**

``` r
print(varpart_ITS$varpart)
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = otu_rclr, X = farm, soil, clim)
    ## 
    ## Explanatory tables:
    ## X1:  farm
    ## X2:  soil
    ## X3:  clim 
    ## 
    ## No. of explanatory tables: 3 
    ## Total variation (SS): 50404 
    ##             Variance: 638.03 
    ## No. of observations: 80 
    ## 
    ## Partition table:
    ##                       Df R.square Adj.R.square Testable
    ## [a+d+f+g] = X1        11  0.37029      0.26842     TRUE
    ## [b+d+e+g] = X2        10  0.72127      0.68088     TRUE
    ## [c+e+f+g] = X3         2  0.38074      0.36465     TRUE
    ## [a+b+d+e+f+g] = X1+X2 21  0.80465      0.73392     TRUE
    ## [a+c+d+e+f+g] = X1+X3 13  0.54452      0.45481     TRUE
    ## [b+c+d+e+f+g] = X2+X3 12  0.81303      0.77954     TRUE
    ## [a+b+c+d+e+f+g] = All 23  0.87072      0.81762     TRUE
    ## Individual fractions                                   
    ## [a] = X1 | X2+X3      11               0.03808     TRUE
    ## [b] = X2 | X1+X3      10               0.36282     TRUE
    ## [c] = X3 | X1+X2       2               0.08371     TRUE
    ## [d]                    0               0.05207    FALSE
    ## [e]                    0               0.10267    FALSE
    ## [f]                    0               0.01496    FALSE
    ## [g]                    0               0.16331    FALSE
    ## [h] = Residuals                        0.18238    FALSE
    ## Controlling 1 table X                                  
    ## [a+d] = X1 | X3       11               0.09015     TRUE
    ## [a+f] = X1 | X2       11               0.05304     TRUE
    ## [b+d] = X2 | X3       10               0.41489     TRUE
    ## [b+e] = X2 | X1       10               0.46549     TRUE
    ## [c+e] = X3 | X1        2               0.18638     TRUE
    ## [c+f] = X3 | X2        2               0.09866     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(varpart_ITS$varpart, cutoff = 0, digits = 1, bg = c("#E1C901", "#DC3220", "#1E88E5"))
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_venn_plot_ITS-1.png)<!-- -->

**TOTAL contribution per variable group**

| Group             | Percentage |
|:------------------|-----------:|
| Farming practices |       35.6 |
| Soil              |       59.3 |
| Climate           |       36.5 |

## 5.3 Cabon cycling functions

**Variance partitionning :**

``` r
varpart_C_fct <- varpart_analysis(otu_C_fct, env_fct)
```

**Variables groups + global significativity :**

``` r
knitr::kable(varpart_C_fct$signifs)
```

| Group             | model_p | model_padj |
|:------------------|--------:|-----------:|
| Global model      |  0.0001 |     0.0004 |
| Farming practices |  0.0207 |     0.0414 |
| Soil              |  0.0001 |     0.0004 |
| Climate           |  0.3988 |     0.3988 |

**varpart() output**

``` r
print(varpart_C_fct$varpart)
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = otu_rclr, X = farm, soil, clim)
    ## 
    ## Explanatory tables:
    ## X1:  farm
    ## X2:  soil
    ## X3:  clim 
    ## 
    ## No. of explanatory tables: 3 
    ## Total variation (SS): 4.25 
    ##             Variance: 0.055922 
    ## No. of observations: 77 
    ## 
    ## Partition table:
    ##                       Df R.square Adj.R.square Testable
    ## [a+d+f+g] = X1        11  0.20826      0.07427     TRUE
    ## [b+d+e+g] = X2        10  0.64372      0.58973     TRUE
    ## [c+e+f+g] = X3         2  0.16066      0.13798     TRUE
    ## [a+b+d+e+f+g] = X1+X2 21  0.73126      0.62866     TRUE
    ## [a+c+d+e+f+g] = X1+X3 13  0.33450      0.19717     TRUE
    ## [b+c+d+e+f+g] = X2+X3 12  0.65330      0.58829     TRUE
    ## [a+b+c+d+e+f+g] = All 23  0.74149      0.62931     TRUE
    ## Individual fractions                                   
    ## [a] = X1 | X2+X3      11               0.04101     TRUE
    ## [b] = X2 | X1+X3      10               0.43214     TRUE
    ## [c] = X3 | X1+X2       2               0.00065     TRUE
    ## [d]                    0               0.01818    FALSE
    ## [e]                    0               0.12225    FALSE
    ## [f]                    0              -0.00209    FALSE
    ## [g]                    0               0.01717    FALSE
    ## [h] = Residuals                        0.37069    FALSE
    ## Controlling 1 table X                                  
    ## [a+d] = X1 | X3       11               0.05919     TRUE
    ## [a+f] = X1 | X2       11               0.03892     TRUE
    ## [b+d] = X2 | X3       10               0.45031     TRUE
    ## [b+e] = X2 | X1       10               0.55439     TRUE
    ## [c+e] = X3 | X1        2               0.12290     TRUE
    ## [c+f] = X3 | X2        2              -0.00144     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(varpart_C_fct$varpart, 
     cutoff = 0, digits = 1, 
     bg = c("#E1C901", "#DC3220", "#1E88E5"))
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_venn_plot_C_cycling-1.png)<!-- -->

**TOTAL contribution per variable group**

| Group             | Percentage |
|:------------------|-----------:|
| Farming practices |       19.9 |
| Soil              |       46.5 |
| Climate           |       13.8 |

## 5.4 Final plot

Saving the necessary data to plot results in the github :

``` r
write.csv(varpart_16S$results_plot, 
          file.path(data_path, "variance_partitionning_16S_results.csv"))
write.csv(varpart_ITS$results_plot, 
          file.path(data_path, "variance_partitionning_ITS_results.csv"))
write.csv(varpart_C_fct$results_plot, 
          file.path(data_path, "variance_partitionning_carbon_cycling_results.csv"))
```

Load the data using this next chunk :

``` r
varpart_16S <- list()
varpart_ITS <- list()
varpart_C_fct <- list()

# dbRDA-associated analyses from this script
varpart_16S$results_plot   <- read.csv(
  file.path(data_path, "variance_partitionning_16S_results.csv"), 
  stringsAsFactors = TRUE)
varpart_ITS$results_plot   <- read.csv(
  file.path(data_path, "variance_partitionning_ITS_results.csv"), 
  stringsAsFactors = TRUE)
varpart_C_fct$results_plot <- read.csv(
  file.path(data_path, "variance_partitionning_carbon_cycling_results.csv"), 
  stringsAsFactors = TRUE)
```

Creating individual plots using the pre-defined function :

``` r
# Making the individual plots
patch_16S_varpart   <- plot_varpart(varpart_16S$results_plot, "(a) 16S ASVs")
patch_ITS_varpart   <- plot_varpart(varpart_ITS$results_plot, "(c) ITS ASVs")
patch_C_fct_varpart <- plot_varpart(varpart_C_fct$results_plot, "(b) Carbon cycling")
```

``` r
# Assembling the plots
patch_16S_varpart + patch_C_fct_varpart + patch_ITS_varpart + guide_area() +
  plot_layout(guides = 'collect') &
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size = 14))
```

![](C:/Users/Jeann/Documents/manuscrit/microbiome-ESSAQ_article/figures/01_fig3-1.png)<!-- -->

``` r
ggsave("fig3.pdf", device = "pdf", path = article_fig_path,
       width = 8, height = 8, units = c("in"))
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
    ##  [1] lubridate_1.9.5   forcats_1.0.1     stringr_1.6.0     dplyr_1.2.1      
    ##  [5] purrr_1.2.2       readr_2.2.0       tidyr_1.3.2       tibble_3.3.1     
    ##  [9] tidyverse_2.0.0   skimr_2.2.2       vegan_2.7-5       permute_0.9-10   
    ## [13] microeco_2.3.0    wesanderson_0.3.7 patchwork_1.3.2   ggpubr_1.0.0     
    ## [17] ggplot2_4.0.3     knitr_1.51        conflicted_1.2.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1    olsrr_0.7.0         farver_2.1.2       
    ##  [4] S7_0.2.2            fastmap_1.2.0       digest_0.6.39      
    ##  [7] timechange_0.4.0    lifecycle_1.0.5     cluster_2.1.8.3    
    ## [10] magrittr_2.0.5      compiler_4.4.1      rlang_1.3.0        
    ## [13] tools_4.4.1         igraph_2.3.3        yaml_2.3.12        
    ## [16] corrplot_0.95       data.table_1.18.6.1 ggsignif_0.6.4     
    ## [19] labeling_0.4.3      plyr_1.8.9          repr_1.1.7         
    ## [22] RColorBrewer_1.1-3  abind_1.4-8         withr_3.0.3        
    ## [25] grid_4.4.1          scales_1.4.0        MASS_7.3-66        
    ## [28] cli_3.6.6           rmarkdown_2.31      ragg_1.5.2         
    ## [31] generics_0.1.4      otel_0.2.0          rstudioapi_0.19.0  
    ## [34] reshape2_1.4.5      tzdb_0.5.0          ape_5.8-1          
    ## [37] cachem_1.1.0        splines_4.4.1       parallel_4.4.1     
    ## [40] BiocManager_1.30.27 base64enc_0.1-6     vctrs_0.7.3        
    ## [43] Matrix_1.7-6        jsonlite_2.0.0      carData_3.0-6      
    ## [46] car_3.1-5           hms_1.1.4           rstatix_1.1.0      
    ## [49] Formula_1.2-6       systemfonts_1.3.2   nortest_1.0-4      
    ## [52] goftest_1.2-3       glue_1.8.1          codetools_0.2-20   
    ## [55] cowplot_1.2.0       stringi_1.8.9       gtable_0.3.6       
    ## [58] pillar_1.11.1       htmltools_0.5.9     R6_2.6.1           
    ## [61] textshaping_1.0.5   evaluate_1.0.5      lattice_0.23-1     
    ## [64] backports_1.5.1     memoise_2.0.1       broom_1.0.13       
    ## [67] renv_1.1.4          Rcpp_1.1.2          gridExtra_2.3.1    
    ## [70] nlme_3.1-170        mgcv_1.9-4          xfun_0.60          
    ## [73] pkgconfig_2.0.3
