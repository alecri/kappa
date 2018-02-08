# kappa

### Description
This repository contains the files to typset my doctoral thesis (or kappa) titled **"Novel methods for dose--response meta-analysis"** (Karolinska Institutet, 2018). It is based on the tex template written by [Andrea Discacciati](https://github.com/anddis/phd-thesis).

### Compiling instruction
To compile the thesis and generate the final PDF (*kappa.pdf*) follows these steps:
- Download and install the last version of [`R`](https://www.r-project.org/) and [`Rstudio`](https://www.rstudio.com/).
- Install the [`packrat`](https://rstudio.github.io/packrat/) package by typing `install.packages("packrat")` either in R or Rstudio.
- Get the [TeX distribution](https://www.latex-project.org/get/#tex-distributions) for your OS.
- Clone or download locally this [repository](https://github.com/alecri/kappa.git) and unzip it.
- Oper Rstudio and create a new project based on an existing folder, i.e. the folder where you unzipped the project.
- Open the “Project Options” from the menu Tools. In the panel Sweave, change the program defaults from Sweave to knitr.
- Type `packrat::restore()` in the Rstudio console.
- Oper the file `kappa.Rnw` in Rstudio and press the button “Compile PDF” (or press cmd+shift+enter on Mac or ctrl+shift+enter on Windows/Linux).

OBS: If you get the following error: `Error in ensurePackageSymlink(source, target) ` run type `unlink("packrat/lib-R", recursive = TRUE)` in the Rstudio console, and then restart your R session.


### Structure
The thesis is written combining [`R`](https://www.r-project.org/) and [LaTex](http://www.latex-project.org/) by generating the tex files with [`knitr`](https://yihui.name/knitr/).

The thesis *Rnw* has been deviding into multiple *Rnw* input files. The root/master file is `kappa.Rnw`.
The main folders are:
- chapters: different *Rnw*s chapters of the kappa
- additional: useful *tex* files, included *bib* library
- style: styles for latex sources
- code: a sample of R scripts for different analyses (some included in the chapters)
- data: data (sets) used in the thesis and in the scripts in the code folder
- packrat: for handling `R` package dependencies and enhancing reproducibility
- figure: pictures and time consuming figures
- cache: figures and results to speed up the compilation (change `cache = TRUE` in the `kappa.Rnw`)

### Author
Alessio Crippa, Karolinska Institutet, Stockholm, Sweden

### Session Info

``` r
R> sessionInfo()
R version 3.4.3 (2017-11-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.3

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] hetmeta_0.1.0      metafor_2.0-0      Matrix_1.2-12      rms_5.1-2         
 [5] SparseM_1.77       Hmisc_4.1-1        Formula_1.2-2      survival_2.41-3   
 [9] lattice_0.20-35    kableExtra_0.7.0   scales_0.5.0.9000  gridExtra_2.3     
[13] cowplot_0.9.2      xtable_1.8-2       Epi_2.24           rworldmap_1.3-6   
[17] sp_1.2-7           lubridate_1.7.1    dosresmeta_2.1.0   mvmeta_0.4.7      
[21] forcats_0.2.0      stringr_1.2.0      dplyr_0.7.4        purrr_0.2.4       
[25] readr_1.1.1        tidyr_0.8.0        tibble_1.4.2       ggplot2_2.2.1.9000
[29] tidyverse_1.2.1    knitr_1.19        

loaded via a namespace (and not attached):
 [1] nlme_3.1-131        cmprsk_2.2-7        RColorBrewer_1.1-2  httr_1.3.1         
 [5] rprojroot_1.3-2     numDeriv_2016.8-1   tools_3.4.3         backports_1.1.2    
 [9] R6_2.2.2            rpart_4.1-12        lazyeval_0.2.1      colorspace_1.3-2   
[13] nnet_7.3-12         mnormt_1.5-5        compiler_3.4.3      cli_1.0.0          
[17] rvest_0.3.2         quantreg_5.34       htmlTable_1.11.2    xml2_1.2.0         
[21] sandwich_2.4-0      checkmate_1.8.5     polspline_1.1.12    mvtnorm_1.0-7      
[25] psych_1.7.8         digest_0.6.15       foreign_0.8-69      rmarkdown_1.8      
[29] base64enc_0.1-3     pkgconfig_2.0.1     htmltools_0.3.6     maps_3.2.0         
[33] htmlwidgets_1.0     rlang_0.1.6         readxl_1.0.0        rstudioapi_0.7     
[37] bindr_0.1           zoo_1.8-1           jsonlite_1.5        acepack_1.4.1      
[41] magrittr_1.5        dotCall64_0.9-5.2   Rcpp_0.12.15        munsell_0.4.3      
[45] stringi_1.1.6       multcomp_1.4-8      yaml_2.1.16         MASS_7.3-48        
[49] plyr_1.8.4          grid_3.4.3          maptools_0.9-2      parallel_3.4.3     
[53] crayon_1.3.4        haven_1.1.1         splines_3.4.3       hms_0.4.1          
[57] pillar_1.1.0        codetools_0.2-15    reshape2_1.4.3      glue_1.2.0         
[61] packrat_0.4.8-1     evaluate_0.10.1     latticeExtra_0.6-28 data.table_1.10.4-3
[65] modelr_0.1.1        spam_2.1-2          MatrixModels_0.4-1  cellranger_1.1.0   
[69] gtable_0.2.0        etm_0.6-2           assertthat_0.2.0    broom_0.4.3        
[73] viridisLite_0.2.0   fields_9.6          bindrcpp_0.2        cluster_2.0.6      
[77] TH.data_1.0-8    
```
