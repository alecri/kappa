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
