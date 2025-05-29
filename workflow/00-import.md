00-import
================
Compiled at 2025-05-29 17:10:33 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "716af71f-d7f8-4fa3-b34f-babf5d3a1051")
```

The purpose of this document is to copy the data files from the raw data
directory to the target directory, so that they can be used in the
workflow. The files copied include CPI scores and Map files.

``` r
library("tidyverse")
```

    ## Warning: package 'ggplot2' was built under R version 4.3.2

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## CPI scores

``` r
INDIR <- "../CPI_scores"

score_files <- list.files(INDIR, full.names = TRUE)
score_files
```

    ## [1] "../CPI_scores/Map.txt"          "../CPI_scores/Scores_dmarA.txt"
    ## [3] "../CPI_scores/Scores_drob.txt"  "../CPI_scores/Scores_dsoxS.txt"
    ## [5] "../CPI_scores/Scores_WT.txt"

``` r
file.copy(score_files, path_target())
```

    ## [1] TRUE TRUE TRUE TRUE TRUE

## Files written

These files have been written to the target directory, `data/00-import`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 5 × 4
    ##   path             type         size modification_time  
    ##   <fs::path>       <fct> <fs::bytes> <dttm>             
    ## 1 Map.txt          file          26K 2025-05-29 17:10:34
    ## 2 Scores_WT.txt    file         106K 2025-05-29 17:10:34
    ## 3 Scores_dmarA.txt file         105K 2025-05-29 17:10:34
    ## 4 Scores_drob.txt  file         104K 2025-05-29 17:10:34
    ## 5 Scores_dsoxS.txt file         105K 2025-05-29 17:10:34

## Session Info

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS 15.4.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.2    
    ##  [5] purrr_1.0.2     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1   
    ##  [9] ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.4        crayon_1.5.2        compiler_4.3.1     
    ##  [4] tidyselect_1.2.0    scales_1.3.0        yaml_2.3.7         
    ##  [7] fastmap_1.1.1       here_1.0.1          R6_2.5.1           
    ## [10] generics_0.1.3      knitr_1.43          munsell_0.5.0      
    ## [13] rprojroot_2.0.3     pillar_1.9.0        tzdb_0.4.0         
    ## [16] rlang_1.1.4         utf8_1.2.4          stringi_1.7.12     
    ## [19] xfun_0.40           fs_1.6.3            timechange_0.2.0   
    ## [22] cli_3.6.3           withr_2.5.2         magrittr_2.0.3     
    ## [25] digest_0.6.33       grid_4.3.1          rstudioapi_0.15.0  
    ## [28] hms_1.1.3           lifecycle_1.0.4     vctrs_0.6.5        
    ## [31] evaluate_0.21       glue_1.7.0          projthis_0.0.0.9025
    ## [34] fansi_1.0.6         colorspace_2.1-0    rmarkdown_2.24     
    ## [37] tools_4.3.1         pkgconfig_2.0.3     htmltools_0.5.6
