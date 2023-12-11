# Handling missing disease information due to death in diseases that need two visits to diagnose
This repository contains the code for the paper entitled: **Handling missing disease information due to death in diseases that need two visits to diagnose** by Le Thi Phuong Thao, Rory Wolfe, Stephane Heritier, and Ronald Geskus, submitted to the Statistics in Medicine journal. 

The author mainly responsible for writing the code and whom readers should approach with questions or bug reports is [Le Thi Phuong Thao](mailto:thao.le@monash.edu).
The repo contains the coding scripts and results from the simulation study described in the paper as follows:

The folder `Simulation` has the following contents:

-   `01_Functions.R`: supporting functions for the simulation
-   `02_RunFuns.R`: functions to run simulations
-   `03_RunSim_Bias_len05.R`, and `03_RunSim_Bias_len10.R`: load the supporting functions and run the simulation for considered scenarios in the manuscript (visit interval = 0.5 and visit interval = 1, respectively). Outputs include the estimate of beta01, shape and scale parameters of Weibul distribution, and number of events for each scenario.
-   `03_RunSim_Bias_Nevent.R`: load the supporting functions, run the simulation and output the number of false positive events in case the false positive rate = 1
-   `04_Summary-Simulation-Results.R`: This R script summaries simulation results, and output the following statistics for each scenario:
    -   Number of events
    -   Estimated bias, and 95% confidence interval
    -   Coverage

The folder `Case study` has the following contents:

-   `DataApplication.R`: This R script does the following task:
    -   Derive the ready to work data set from the ASPREE data. Note that data can be obtained via the AMS system of ASPREE (<https://ams.aspree.org/public/>)
    -   Fit multiple regression models using ASPREE data
    -   Output results

The folder `Interdata` contains the corresponding R work-spaces with the simulation results, summary of simulation results, and case study results:

-   `2023-02-17_simresults-All-Length05.Rdata`: Output of `03_RunSim_Bias_len05.R`

-   `2023-02-17_simresults-All-Length05.Rdata`: Output of `03_RunSim_Bias_len10.R`

-   `2023-02-20_simresults-Fasel-Positive-Rate.Rdata`: Output of `03_RunSim_Bias_Nevent.R`

-   `2023-02-20_Bias-Results.Rdata`: Output of `04_Summary-Simulation-Results.R`

-   `2023-02-20_Nevents-Results.Rdata`: Output of `04_Summary-Simulation-Results.R`

-   `2023-02-20_Simulation-Senarios.Rdata`: Output of `04_Summary-Simulation-Results.R`

-   `2023-01-18_aspree-idm.Rdata`: Output of `DataApplication.R`

The folder `Results` contains:

-   `TablesFigures.Rmd`: an R markdown file to generate figures and tables in the manuscript
-   `TablesFigures.html`: the corresponding html output.

The results have been produced in R under the following specification:

```         
R version 4.0.5 (2021-03-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/local/intel/2018u3/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so

locale:
 [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8    
 [5] LC_MONETARY=en_AU.UTF-8    LC_MESSAGES=en_AU.UTF-8   
 [7] LC_PAPER=en_AU.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggh4x_0.2.6        RColorBrewer_1.1-2 wesanderson_0.3.7  here_1.0.1        
 [5] gtsummary_1.7.2    broom_1.0.5        glue_1.6.2         gt_0.10.0         
 [9] forcats_1.0.0      stringr_1.4.0      readr_1.4.0        tibble_3.2.1      
[13] ggplot2_3.4.4      tidyverse_1.3.1    icenReg_2.0.15     coda_0.19-4       
[17] Rcpp_1.0.7         SurvRegCensCov_1.4 mstate_0.3.1       furrr_0.2.3       
[21] future_1.21.0      SmoothHazard_1.4.1 prodlim_2019.11.13 purrr_1.0.2       
[25] tidyr_1.1.4        dplyr_1.1.4        survival_3.2-13   

loaded via a namespace (and not attached):
 [1] httr_1.4.2           jsonlite_1.7.2       splines_4.0.5       
 [4] foreach_1.5.1        modelr_0.1.8         assertthat_0.2.1    
 [7] cellranger_1.1.0     globals_0.14.0       numDeriv_2016.8-1.1 
[10] pillar_1.9.0         backports_1.2.1      lattice_0.20-44     
[13] digest_0.6.27        rvest_1.0.0          colorspace_2.0-2    
[16] htmltools_0.5.7      Matrix_1.3-4         pkgconfig_2.0.3     
[19] listenv_0.8.0        haven_2.5.4          scales_1.3.0        
[22] lava_1.6.9           generics_0.1.0       ellipsis_0.3.2      
[25] withr_2.5.2          cli_3.6.1            magrittr_2.0.3      
[28] crayon_1.4.1         readxl_1.3.1         fs_1.6.3            
[31] fansi_0.5.0          parallelly_1.26.1    broom.helpers_1.14.0
[34] xml2_1.3.6           tools_4.0.5          data.table_1.14.0   
[37] hms_1.1.0            lifecycle_1.0.4      munsell_0.5.0       
[40] reprex_2.0.0         compiler_4.0.5       rlang_1.1.2         
[43] grid_4.0.5           iterators_1.0.13     rstudioapi_0.13     
[46] gtable_0.3.0         codetools_0.2-18     DBI_1.1.1           
[49] R6_2.5.0             lubridate_1.7.10     knitr_1.45          
[52] fastmap_1.1.0        utf8_1.2.1           rprojroot_2.0.2     
[55] stringi_1.7.3        parallel_4.0.5       vctrs_0.6.5         
[58] dbplyr_2.1.1         tidyselect_1.2.0     xfun_0.41           
```
